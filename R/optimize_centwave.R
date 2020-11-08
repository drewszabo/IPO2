#' @import data.table
#'
optimize_centwave <- function(
  raw_data = NULL,
  parameter_list = suggest_centwave_params(),
  bpparam = BiocParallel::bpparam(),
  log_file = NULL,
  plot_dir = NULL
) {

  # check parameters
  check_centwave_params(parameter_list)

  # create plot directory
  pd <- paste0(plot_dir, "/centwave_contours/")
  if (!is.null(plot_dir)) {
    if (!dir.exists(plot_dir)) {
      stop("The plot directory, ", plot_dir, "does not exist")
    } else if (!dir.exists(pd)) {
      dir.create(pd)
    }
  }

  # set up iteration loop
  history <- list()
  iteration <- 1
  while(iteration <= 50) {

    message(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      " : Starting iteration ",
      iteration
    )

    # parse parameters
    parameters <- parse_parameters(parameter_list)

    # design experiments
    design <- design_experiments(parameters)

    # generate centWave parameters
    cwp <-
      purrr::pmap(
        cbind(design, t(parameters$constant), stringsAsFactors = FALSE),
        make_cwp,
        roiList = parameters$lists$roiList,
        roiScales = parameters$lists$roiScales
      )

    # run xcms for each experiment
    xcmsnexp <-
      BiocParallel::bplapply(
        seq_along(cwp),
        function(x) {
          is_try_error <- TRUE
          j <- 1
          while(is_try_error & j <= 5) {
            xcmsnexp_trial <- try(
              xcms::findChromPeaks(
                raw_data,
                param = cwp[[x]],
                BPPARAM = BiocParallel::SerialParam())
            )
            is_try_error <- test_try_error(xcmsnexp_trial)
            if (!is.null(log_file)) {
              sink(log_file, append = TRUE, type = "output")
              cat("Iteration:", sprintf("%02d", iteration),
                  "   CWP:", sprintf("%02d", x),
                  "   Attempt:", j,
                  "   Try error:", is_try_error, "\n")
            }
            j <- j + 1
          }
          xcmsnexp_trial
        }
      )

    # score experiment
    score <-
      rbindlist(
        BiocParallel::bplapply(
          xcmsnexp,
          score_peaks
        )
      )

    # calculate model
    model <- create_model(design, score)

    # find best parameters
    maximum <- get_maximum(design, model)

    # make plots
    if (!is.null(plot_dir)) {
      plot_name <- paste0(pd, "centwave_contour_", iteration, ".png")
      plot_contours(design, model, maximum, plot_name)
    }

    # score best parameters
    matched_row <- design[as.list(maximum), on = names(maximum), which = TRUE]

    if(!is.na(matched_row)) {
      max_score <- score[matched_row, ]
      max_cwp <- cwp[matched_row]
    } else {
      max_cwp <- purrr::pmap(
        c(maximum, parameters$constant),
        make_cwp,
        roiList = parameters$lists$roiList,
        roiScales = parameters$lists$roiScales
      )
      max_xcmsnexp <- xcms::findChromPeaks(raw_data, param = max_cwp[[1]])
      max_score <- score_peaks(max_xcmsnexp)
    }

    # assign to history
    history[[iteration]] <-
      cbind(
        t(unlist(parameters$to_optimize)),
        cwp = max_cwp,
        t(maximum),
        max_score
      )

    # check for improvement
    if (iteration == 1) {
      better <- TRUE
    } else {
      better <- history[[iterator]][["score"]] > history[[iterator - 1]][["score"]]
    }

    # adjust intervals
    if (better) {
      new_parameters <- pick_parameters(parameters$to_optimize, maximum)
      parameter_list[names(new_parameters)] <- new_parameters
    } else {
      break
    }

    # increase counter
    iteration <- iteration + 1

  }

  # output results
  history <- rbindlist(history)
  best_params <- history[score == max(score), cwp]
  cat(best_params)
  list(history = history, best_params = best_params)

}


get_centwave_defaults <- function() {
  list(
    ppm                = 25,
    min_peakwidth      = 20,
    max_peakwidth      = 50,
    snthresh           = 10,
    prefilter_k        = 3,
    prefilter_int      = 100,
    mzdiff             = -0.001,
    noise              = 0,
    mzCenterFun        = "wMean",
    integrate          = 1L,
    fitgauss           = FALSE,
    verboseColumns     = FALSE,
    roiList            = list(),
    firstBaselineCheck = TRUE,
    roiScales          = numeric()
  )
}


suggest_centwave_params <- function() {
  list(
    ppm                = c(17, 32),
    min_peakwidth      = c(12, 28),
    max_peakwidth      = c(35, 65),
    snthresh           = 100,  # originally 10, increase to increase speed
    prefilter_k        = 3,
    prefilter_int      = 1000,  # originally 100
    mzdiff             = c(-0.001, 0.01),
    noise              = 0  # originally 0, increase to increase speed
  )
}


param_types <- function() {

  assign(
    "quant",
    c(
      "ppm",
      "min_peakwidth",
      "max_peakwidth",
      "snthresh",
      "prefilter_k",
      "prefilter_int",
      "mzdiff",
      "noise"
    ),
    envir = parent.frame()
  )

  assign(
    "qual",
    c(
      "mzCenterFun",
      "integrate",
      "fitgauss",
      "verboseColumns",
      "firstBaselineCheck"
    ),
    envir = parent.frame()
  )

  assign(
    "lists",
    c(
      "roiList",
      "roiScales"
    ),
    envir = parent.frame()
  )

}


check_centwave_params <- function(parameter_list) {

  param_types()

  # check parameter names
  nms <- names(parameter_list)
  bad_names <- nms[nms %nin% c(quant, qual, lists)]
  if (length(bad_names) != 0) {
    stop("These parameters are not used by centWave: ",
         paste(bad_names, collapse = ", "))
  }

  # check if too many parameters provided
  lengths <- lapply(parameter_list[quant], length)
  bad_lengths <- names(lengths[lengths > 2])
  if (any(lengths > 2)) {
    stop("Too many parameters specified for: ",
         paste(bad_lengths, collapse = ", "))
  }

  # check for parameters to optimize
  if (!(any(lengths == 2))) {
    stop("No parameters specified for optimization")
  }

  # check ranges
  for (nm in nms){
    if (nm %in% c("ppm")) {
      if (min(parameter_list[[nm]]) <= 0) {
        stop(paste("The parameter", nm, "must be greater than 0"))
      }
    } else if (nm %in% c("min_peakwidth", "max_peakwidth")) {
      if (min(parameter_list[[nm]]) <= 3) {
        stop(paste("The parameter", nm, "must be greater than 3"))
      }
    } else if (nm %in% c("snthresh", "prefilter_k", "prefilter_int", "noise")) {
      if (min(parameter_list[[nm]]) < 0) {
        stop(paste("The parameter", nm, "must not be negative"))
      }
    }
  }
}


parse_parameters <- function(parameters) {
  parameter_list <- get_centwave_defaults()
  parameter_list[names(parameters)] <- parameters
  param_types()
  params <- list()
  params$to_optimize <- parameter_list[quant][lapply(parameter_list[quant], length) == 2]
  params$constant <- parameter_list[c(quant, qual)][lapply(parameter_list[c(quant, qual)], length) != 2]
  params$lists <- parameter_list[lists]
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


make_cwp <- function(
  min_peakwidth = 20,
  max_peakwidth = 50,
  prefilter_k = 3,
  prefilter_int = 100,
  ...
) {
  xcms::CentWaveParam(
    peakwidth = c(min_peakwidth, max_peakwidth),
    prefilter = c(prefilter_k, prefilter_int),
    ...)
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
  number_columns <- c(1, 1, 1, 2, 2, 3, 3, 4)
  number_rows <- number_pairs / number_columns

  png(
    plot_name,
    width = 7.5,
    height = 10,
    units = "in",
    res = 300
  )

  params <- ncol(design)
  cols <- number_columns[which(number_params == params)]
  rows <- number_rows[which(number_params == params)]
  par(mfrow = c(rows, cols))

  # plot
  contour(model, form = form, image = TRUE, at = maximum)
  dev.off()

}


pick_parameters <- function(parameters, maximum) {

  radius <- purrr::map_dbl(parameters, ~diff(.x) / 2)
  params <- list()

  for (nm in names(parameters)) {

    if (min(parameters[[nm]]) == maximum[nm]) {
      upper <- maximum[nm] + 1.2 * radius[nm]
      lower <- maximum[nm] - 1.2 * radius[nm]
    } else if (max(parameters[[nm]]) == maximum[nm]) {
      lower <- maximum[nm] - 1.2 * radius[nm]
      upper <- maximum[nm] + 1.2 * radius[nm]
    } else {
      lower <- maximum[nm] - radius[nm] * 0.8
      upper <- maximum[nm] + radius[nm] * 0.8
    }

    if (nm %in% c("ppm")) {
      lower <- max(1, lower)
    } else if (nm %in% c("min_peakwidth", "max_peakwidth")) {
      lower <- max(3, lower)
    } else if (nm %in% c("snthresh", "noise", "prefilter_k", "prefilter_int")) {
      lower <- max(0, lower)
    }

    params[[nm]] <- unique(c(lower, upper))  # stop optimizing if converge

  }

  params

}
