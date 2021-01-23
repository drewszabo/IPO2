"%nin%" <- function(x, table) {
  match(x, table, nomatch = 0L) == 0L
}

prepare_out_dir <- function(out_dir = NULL, fun_name) {

  if (is.null(out_dir)) out_dir <- getwd()
  if(!dir.exists(out_dir)) {
    cat("Creating", out_dir, "\n")
    dir.create(out_dir, recursive = TRUE)
  }

  dir <- paste0(out_dir, "/", fun_name, format(Sys.time(), "_%Y-%m-%d-%H-%M-%S"))

  dir.create(dir)
  assign("dir", paste0(dir, "/"), envir = parent.frame())
  assign("log_file", paste0(dir, "/log_file.txt"), envir = parent.frame())

}

parse_parameters <- function(
  parameter_list,
  method = c("centwave", "align_group")
) {

  default_params <- switch(
    method,
    "centwave" = {
      param_types_centwave()
      get_centwave_defaults()
    },
    "align_group" = {
      param_types_align_group()
      c(get_obiwarp_defaults(), get_density_defaults())
    }
  )

  default_params[names(parameter_list)] <- parameter_list

  params <- list()
  params$to_optimize <- default_params

  params <- list()
  params$to_optimize <- Filter(function(x) length(x) == 2, default_params)
  params$constant <- Filter(function(x) length(x) != 2, default_params[c(quant, qual)])
  params$lists <- default_params[lists]

  params
}


generate_ccd <- function(parameters_to_optimize) {

  lower_bound <- sapply(parameters_to_optimize, "[[", 1)
  upper_bound <- sapply(parameters_to_optimize, "[[", 2)
  center <- (upper_bound - lower_bound) / 2

  x <- paste0(
    "x",
    1:length(parameters_to_optimize),
    " ~ (",
    c(names(parameters_to_optimize)),
    " - ",
    (lower_bound + center),
    ") / ",
    center
  )

  formulas <- lapply(x, as.formula)

  rsm::ccd(
    length(parameters_to_optimize),  # number of variables
    n0 = 1,  # number of center points
    alpha = "face",  # position of the ‘star’ points
    randomize = FALSE,
    inscribed = TRUE,  # axis points are at +/- 1 and the cube points are at interior positions
    coding = formulas  # list of coding formulas for the design variables
  )

}


design_experiments <- function(parameters) {

  if(length(parameters$to_optimize) > 1) {
    design <- generate_ccd(parameters$to_optimize)
    design <- rsm::decode.data(design)
    design <- subset(design, select = -c(run.order, std.order, Block))
  } else {
    design <- seq(min(parameters$to_optimize[[1]]),
                  max(parameters$to_optimize[[1]]),
                  length.out = 9)
    design <- data.frame(design)
    colnames(design) <- names(parameters$to_optimize[1])
  }

  data.table(design)

}


create_model <- function(design, score) {

  params <- paste0(colnames(design), collapse = ", ")
  design <- cbind(design, score = score$score)

  if(ncol(design) > 1) {
    formula <- as.formula(paste0("score ~ SO(", params, ")"))
    model <- rsm::rsm(formula, data = design)
  } else {
    formula <- as.formula(paste("score ~", params, "+", params, "^ 2"))
    model <- lm(formula, data = design)
  }

  model

}


get_maximum <- function(design, model) {

  number_params <- ncol(design)
  steps <- ceiling(1E6 ^ (1/number_params))

  param_values <- list()
  for (nm in names(design)) {
    minimum <- min(design[[nm]])
    maximum <- max(design[[nm]])
    param_values[[nm]] <- seq(minimum, maximum, length.out = steps)
  }

  search_grid <- do.call(CJ, param_values)
  max_value <- predict(model, search_grid)
  unlist(search_grid[max_value == max(max_value), ][1, ])

}


plot_contours <- function(design, model, maximum, plot_name) {

  # make formula
  form <- as.formula(
    paste("~ ", paste(colnames(design), collapse = " + "))
  )

  # set plotting area
  number_params <- 1:8
  number_pairs <- number_params * (number_params - 1) / 2
  number_rows <- c(1, 1, 1, 2, 2, 3, 3, 4)
  number_columns <- number_pairs / number_rows

  params <- ncol(design)
  cols <- number_columns[which(number_params == params)]
  rows <- number_rows[which(number_params == params)]

  png(
    plot_name,
    width = cols * 4,
    height = rows * 4,
    units = "in",
    res = 300
  )

  par(mfrow = c(rows, cols))

  # plot
  contour(model, form = form, image = TRUE, at = maximum)
  dev.off()

}


pick_parameters <- function(parameters, maximum) {

  width <- purrr::map_dbl(parameters, diff)
  delta <- purrr::map2(parameters, as.list(maximum), ~abs(.x - .y))
  params <- list()

  boundaries <- list(
    ppm = c(0, Inf, 2),
    min_peakwidth = c(0, Inf, 2),
    max_peakwidth = c(0, Inf, 2),
    snthresh = c(0, Inf, 0),
    noise = c(0, Inf, 0),
    prefilter_k = c(0, Inf, 0),
    prefilter_int = c(0, Inf, 0),
    mzdiff = c(-Inf, Inf, 4),
    binSize_O = c(0, Inf, 3),
    response = c(1, 100, 0),
    gapInit = c(0, Inf, 2),
    gapExtend = c(0, Inf, 2),
    factorDiag = c(-Inf, Inf, 0),
    factorGap = c(-Inf, Inf, 0),
    initPenalty = c(0, Inf, 2),
    bw = c(1, Inf, 2),
    minFraction = c(0, 1, 2),
    minSamples = c(0, Inf, 0),
    binSize_D = c(0.001, Inf, 3),
    maxFeatures = c(0, Inf, 2)
  )

  for (nm in names(parameters)) {

    if (sum(delta[[nm]] == 0) > 0) {
      width[[nm]] <- width[[nm]] * 1.4
    } else {
      width[[nm]] <- width[[nm]] * 0.8
    }

    upper <- maximum[[nm]] + 0.5 * width[[nm]]
    lower <- maximum[[nm]] - 0.5 * width[[nm]]

    lower <- max(boundaries[[nm]][[1]], lower)
    upper <- min(boundaries[[nm]][[2]], upper)

    lower <- round(lower, digits = boundaries[[nm]][[3]])
    upper <- round(upper, digits = boundaries[[nm]][[3]])

    params[[nm]] <- unique(c(lower, upper))  # stop optimizing if converge

  }

  # prevent overlapping peakwidth optimization
  if (sum(
    is.null(params[["min_peakwidth"]]), is.null(params[["max_peakwidth"]])) == 0
  ){
    if (max(params[["min_peakwidth"]] > min(params[["max_peakwidth"]]))) {
      params[["max_peakwidth"]][[1]] <- params[["min_peakwidth"]][[2]]
    }
  }

  params

}


retry_parallel <- function(fun) {
  redo <- TRUE
  trial <- 1
  redo_list <- list()
  while(redo & trial <= 5) {
    out <-
      BiocParallel::bptry(
        eval(fun, env = rlang::env(rlang::caller_env(), redo_list = redo_list))
      )
    errs <- sum(BiocParallel::bpok(out) == FALSE)
    redo <- errs > 0
    if (errs > 0) {
      ids <- which(BiocParallel::bpok(out) == FALSE)
    } else {
      ids <- vector()
    }
    cat(
      "     Trial:", trial,
      "     Errors:", sprintf("%3d", errs),
      "     IDs:", ids,
      "\n"
    )
    if (redo) {
      redo_list <- out
      trial <- trial + 1
    }
  }
  out
}
