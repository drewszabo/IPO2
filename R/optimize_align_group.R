#' Find peak alignment and grouping parameters
#'
#' `optimize_align_group` will identify parameters for `xcms::adjustRtime` using
#' the ObiWarp algorithm and for `xcms::groupChromPeaks` using peak densities.
#'
#' @param xcmsnexp An XCMSnExp object
#' @param parameter_list List of parameters to set or to optimize
#' @param out_dir Path for output directory. Defaults to working directory.
#'
#' @return A list of length 3 containing:
#' \describe{
#'   \item{history}{a data.table object recording the optimization history}
#'   \item{best_obi}{an ObiWarpParam object of the optimal parameters}
#'   \item{best_density}{an PeakDensityParam object of the optimal parameters}
#' }
#' @export
#'
#' @examples
#'
optimize_align_group <- function(
  xcmsnexp = NULL,
  parameter_list = suggest_align_group_params(),
  out_dir = NULL
){

  # prepare output
  prepare_out_dir(out_dir = out_dir, fun_name = "align-group")
  sink(log_file, append = TRUE, type = "output", split = TRUE)
  on.exit(sink(), add = TRUE, after = TRUE)

  # output
  old_w <- getOption("width")
  options(width = 1000)
  on.exit(options(width = old_w), add = TRUE, after = TRUE)

  max_print <- getOption("max.print")
  options(max.print = 99999)
  on.exit(options(max.print = max_print), add = TRUE, after = TRUE)

  # check xcmsnexp
  if (nrow(xcmsnexp) <= 1) {
    stop("Not enough files for alignment")
  }

  # check parameters
  check_align_group_params(parameter_list)
  obi_params <- names(get_obiwarp_defaults())
  density_params <- names(get_density_defaults())

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

    # generate density parameters
    group_params <- params[, colnames(params) %in% density_params, with = FALSE]
    setnames(group_params, "binSize_D", "binSize")
    density <-
      purrr::pmap(
        group_params,
        xcms::PeakDensityParam
      )

    # run alignment
    score <- rbindlist(
      retry_parallel(
        rlang::expr(
          BiocParallel::bpmapply(
            function(x, y) {
              BiocParallel::register(BiocParallel::SerialParam())
              xcms::adjustRtime(
                xcmsnexp,
                param = x
              ) %>%
                xcms::groupChromPeaks(
                  param = y
                ) %>%
                score_align_group()
            },
            x = obi,
            y = density,
            BPREDO = redo_list,
            SIMPLIFY = FALSE
          )
        )
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
    plot_name <-
      paste0(dir, "/contour_align-group_", sprintf("%02d", iteration), ".png")
    plot_contours(design, model, maximum, plot_name)

    # score best parameters
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

    max_score <-
      xcms::adjustRtime(xcmsnexp, param = max_obi[[1]]) %>%
      xcms::groupChromPeaks(param = max_density[[1]]) %>%
      score_align_group()

    # assign to history
    history[[iteration]] <-
      c(
        unlist(parameters$to_optimize),
        obi = list(max_obi),
        density = list(max_density),
        maximum,
        unlist(max_score)
      )

    hx <- rbindlist(history)[, !c("obi", "density")]
    end_cols <- c("rcs", "good", "bad", "gs")
    setcolorder(hx, order(setdiff(names(hx), end_cols)))
    hx[, c("rcs_adj", "gs_adj") := lapply(.SD, scales::rescale, to = c(0, 1)), .SDcols = c("rcs", "gs")
    ][
      , score := rcs_adj + gs_adj
    ]

    cat("\n")
    capture.output(hx, file = log_file, append = TRUE)
    cat("\n\n")

    # check for improvement
    better <- hx[[iteration, "score"]] == max(hx[["score"]])

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
  history[, c("rcs_adj", "gs_adj") := lapply(.SD, scales::rescale, to = c(0, 1)), .SDcols = c("rcs", "gs")
  ][
    , score := rcs_adj + gs_adj
  ]
  setcolorder(history, order(setdiff(names(history), c(end_cols, obi, density))))
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
