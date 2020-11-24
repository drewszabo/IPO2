

optimize_align_group <- function(
  xcmsnexp = NULL,
  parameter_list = suggest_align_group_params(),
  log_file = NULL,
  plot_dir = NULL
){

  # check xcmsnexp
  if (nrow(xcmsnexp) <= 1) {
    stop("Not enough files for alignment")
  }

  # check parameters
  check_align_group_params(parameter_list)
  obi_params <- names(get_obiwarp_defaults())
  density_params <- names(get_density_defaults())

  # log file
  if (!is.null(log_file)) {
    log_file <- check_log_file(log_file)
    sink(log_file, append = TRUE, type = "output")
    on.exit(sink(), add = TRUE, after = TRUE)
  }

  # plot file
  if (!is.null(plot_dir)) {
    pd <- check_plot_dir(plot_dir, "align_group")
  }

  # identify center sample
  dt <- data.table(xcms::chromPeaks(xcmsnexp), keep.rownames = TRUE)
  idx_sample <- dt[, median(intb, na.rm = TRUE), by = sample][, which.max(V1)]
  parameter_list$centerSample <- idx_sample
  if (is.null(parameter_list$sampleGroups)) {
    parameter_list$sampleGroups <- rep(1, length(unique(unlist(dt$sample))))
  }

  # set up iteration loop
  history <- list()
  iteration <- 1
  while(iteration <= 50) {

    cat(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      ": Starting iteration ",
      iteration,
      "\n"
    )

    # parse parameters
    parameters <- parse_parameters(parameter_list, "align_group")

    # design experiments
    design <- design_experiments(parameters)
    params <- cbind(design, t(parameters$constant), t(parameters$lists))

    # generate obiwarp parameters
    align_params <- params[, colnames(params) %in% obi_params, with = FALSE]
    setnames(align_params, "binSize_O", "binSize")
    obi <-
      purrr::pmap(
        align_params,
        xcms::ObiwarpParam
      )

    # run alignment
    redo <- TRUE
    trial <- 1
    redo_list <- list()
    while(redo & trial <= 5) {

      aligned <-
        BiocParallel::bptry(
          BiocParallel::bplapply(
            obi,
            function(x) {
              register(BiocParallel::SerialParam())
              xcms::adjustRtime(
                xcmsnexp,
                param = x
              )
            },
            BPREDO = redo_list
          )
        )

      errs <- sum(bpok(aligned) == FALSE)
      cat(
        "Align Iteration:", sprintf("%02d", iteration),
        "     trial:", trial,
        "     Errors:", errs,
        "\n"
      )
      redo <- errs > 0
      if (redo) {
        redo_list <- aligned
        trial <- trial + 1

      }
    }

    # generate density parameters
    group_params <- params[, colnames(params) %in% density_params, with = FALSE]
    setnames(group_params, "binSize_D", "binSize")
    density <-
      purrr::pmap(
        group_params,
        xcms::PeakDensityParam
      )

    # run grouping
    redo <- TRUE
    trial <- 1
    redo_list <- list()
    while(redo & trial <= 5) {

      grouped <-
        BiocParallel::bptry(
          BiocParallel::bpmapply(
            function(x, y) {
              register(BiocParallel::SerialParam())
              xcms::groupChromPeaks(
                x,
                param = y
              )
            },
            x = aligned,
            y = density,
            BPREDO = redo_list
          )
        )

      errs <- sum(bpok(grouped) == FALSE)
      cat(
        "Grouping Iteration:", sprintf("%02d", iteration),
        "     trial:", trial,
        "     Errors:", errs,
        "\n"
      )
      redo <- errs > 0
      if (redo) {
        redo_list <- grouped
        trial <- trial + 1

      }
    }

    # score
    score <-
      rbindlist(
        BiocParallel::bplapply(
          grouped,
          score_align_group
        )
      )

    score[, c("rcs", "gs") := lapply(.SD, scales::rescale, to = c(0, 1)), .SDcols = c("rcs", "gs")
    ][
      , score := rcs + gs
    ]

    # calculate model
    model <- create_model(design, score)

    # find maximum
    maximum <- get_maximum(design, model)
    max_params <- c(maximum, parameters$constant, parameters$lists)

    # make plots
    if (!is.null(plot_dir)) {
      plot_name <-
        paste0(pd, "align_group_contour_", sprintf("%02d", iteration), ".png")
      plot_contours(design, model, maximum, plot_name)
    }

    # score best parameters
    matched_row <- design[as.list(maximum), on = names(maximum), which = TRUE]

    if(!is.na(matched_row)) {
      max_score <- score[matched_row, ]
      max_obi <- obi[matched_row]
      max_density <- density[matched_row]
    } else {
      align_params <- max_params[names(max_params) %in% obi_params]
      names(align_params)[names(align_params) == "binSize_O"] <- "binSize"
      max_obi <- purrr::pmap(
        data.table(t(align_params)),
        xcms::ObiwarpParam
      )

      group_params <- max_params[names(max_params) %in% density_params]
      names(group_params)[names(group_params) == "binSize_D"] <- "binSize"
      max_density <- purrr::pmap(
        data.table(t(group_params)),
        xcms::PeakDensityParam
      )

      max_xcmsnexp <- xcms::adjustRtime(xcmsnexp, param = max_obi[[1]])
      max_xcmsnexp <- xcms::groupChromPeaks(xcmsnexp, param = max_density[[1]])

      max_score <- score_align_group(max_xcmsnexp)
      max_score <-
        (
          scales::rescale(c(max_score$rcs, score$rcs), to = c(0, 1)) +
            scales::rescale(c(max_score$gs, score$gs), to = c(0, 1))
        )[[1]]

    }

    # assign to history
    history[[iteration]] <-
      c(
        unlist(parameters$to_optimize),
        obi = list(max_obi),
        density = list(max_density),
        maximum,
        score = max_score
      )

    # check for improvement
    if (iteration == 1) {
      better <- TRUE
    } else {
      better <- history[[iteration]][["score"]] > history[[iteration - 1]][["score"]]
    }

    # adjust intervals
    if (better) {
      new_parameters <- pick_parameters(parameters$to_optimize, maximum)
      parameter_list[names(new_parameters)] <- new_parameters
    } else {
      break
    }

    # increment counter
    iteration <- iteration + 1

  }

  # output results
  history <- rbindlist(history)
  best_obi <- history[score == max(score), obi][[1]]
  best_density <- history[score == max(score), density][[1]]
  list(history = history, best_obi = best_obi, best_density = best_density)

}

get_obiwarp_defaults <- function() {
  list(
    binSize_O          = 1,
    centerSample       = integer(),
    response           = 1L, # 0-100; 0 linear on ends, 100 uses all bijective anchors
    distFun            = "cor_opt",
    gapInit            = numeric(), # penalty for gap opening
    gapExtend          = numeric(), # penalty for gap enlargement
    factorDiag         = 2, # local weight for diagonal moves in alignment
    factorGap          = 1, # local wight for gap moves in alignment
    localAlignment     = FALSE, # should a local alignment be performed?
    initPenalty        = 0, # penalty for initiating alignment (local alignment only)
    subset             = integer(), # integer indices of samples for alignment
    subsetAdjust       = "average" # method by which non-subset samples are adjusted
  )
}

get_density_defaults <- function() {
  list(
    sampleGroups       = numeric(), # vector of same length as samples defining groups
    bw                 = 30, # standard deviation of the smooth kernel to be used
    minFraction        = 0.5, # minimum fraction of samples in at least one sample group
    minSamples         = 1, # minimum number of samples in at least one group
    binSize_D          = 0.25, # size of the overlapping slices in m/z dimension
    maxFeatures        = 50 # maximum number of peak groups to be identified in a single slice
  )
}

suggest_align_group_params <- function() {
  list(
    binSize_O          = c(0.7, 1),
    response           = c(1, 10),
    gapInit            = c(0, 0.4),
    gapExtend          = c(2.1, 2.7),
    bw                 = c(22, 38),
    minFraction        = c(0.3, 0.7),
    binSize_D          = c(0.015, 0.035)
  )
}

param_types_align_group <- function() {

  assign(
    "quant",
    c(
      "binSize_O",
      "response",
      "gapInit",
      "gapExtend",
      "factorDiag",
      "factorGap",
      "initPenalty",
      "bw",
      "minFraction",
      "minSamples",
      "binSize_D",
      "maxFeatures"
    ),
    envir = parent.frame()
  )

  assign(
    "qual",
    c(
      "centerSample",
      "distFun",
      "localAlignment",
      "subsetAdjust"
    ),
    envir = parent.frame()
  )

  assign(
    "lists",
    c(
      "subset",
      "sampleGroups"
    ),
    envir = parent.frame()
  )

}

check_align_group_params <- function(parameter_list) {

  param_types_align_group()

  # check parameter names
  nms <- names(parameter_list)
  bad_names <- nms[nms %nin% c(quant, qual, lists)]
  if (length(bad_names) != 0) {
    stop("These parameters are not used: ",
         paste(bad_names, collapse = ", "))
  }

  # check if too many parameters provided
  lengths <- lapply(parameter_list[nms %in% quant], length)
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
  negative <- lapply(parameter_list[nms %in% quant], function(x) min(x) < 0)
  if (any(unlist(negative))) {
    stop("These parameters can not be negative: ",
         paste(parameter_list[unlist(negative)], collapse = ", "))
  }

}
