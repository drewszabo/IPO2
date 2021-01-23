#' Find centWave parameters
#'
#' `optimize_centwave` will identify parameters for `xcms::findChromPeaks` by
#' optimizing the identification of isotope pairs. Note that currently the
#' isotope identification algorithm relies on high resolution m/*z* measurements.
#'
#' @param raw_data MSnExp object containing the raw data
#' @param parameter_list List of parameters to set or to optimize
#' @param out_dir Path for output directory. Defaults to working directory.
#'
#' @return A list of length 2 containing:
#' \describe{
#'   \item{history}{a data.table object recording the optimization history}
#'   \item{best_cwp}{a CentWaveParam object of the optimal parameters}
#' }
#' @export
#'
#' @examples
#'
#'
optimize_centwave <- function(
  raw_data = NULL,
  parameter_list = suggest_centwave_params(),
  out_dir = NULL
) {

  # prepare output
  prepare_out_dir(out_dir = out_dir, fun_name = "centwave")
  sink(log_file, append = TRUE, type = "output", split = TRUE)
  on.exit(sink(), add = TRUE, after = TRUE)

  # check parameters
  check_centwave_params(parameter_list)

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
    parameters <- parse_parameters(parameter_list, "centwave")

    # design experiments
    design <- design_experiments(parameters)

    # generate centWave parameters
    cwp <-
      purrr::pmap(
        cbind(design, t(parameters$constant), t(parameters$lists)),
        make_cwp
      )

    # run xcms for each experiment
    score <- rbindlist(
      retry_parallel(
        rlang::expr(
          BiocParallel::bplapply(
            cwp,
            function(x) {
              score_peaks(
                xcms::findChromPeaks(
                  raw_data,
                  param = x,
                  BPPARAM = BiocParallel::SerialParam()
                )
              )
            },
            BPREDO = redo_list
          )
        )
      )
    )

    # calculate model
    model <- create_model(design, score)

    # find best parameters
    maximum <- get_maximum(design, model)

    # make plots
    plot_name <-
      paste0(dir, "/contour_centwave_", sprintf("%02d", iteration), ".png")
    plot_contours(design, model, maximum, plot_name)

    # score best parameters
    max_cwp <- purrr::pmap(
      data.table(t(maximum), t(parameters$constant), t(parameters$lists)),
      make_cwp
    )
    max_xcmsnexp <- xcms::findChromPeaks(raw_data, param = max_cwp[[1]], BPPARAM = BiocParallel::SerialParam())
    max_score <- score_peaks(max_xcmsnexp)

    # assign to history
    history[[iteration]] <-
      c(
        unlist(parameters$to_optimize),
        cwp = list(max_cwp),
        maximum,
        max_score
      )

    hx <- rbindlist(history)[, !c("cwp")]
    end_cols <- c("total", "isotopes", "indeterminate", "score")
    setcolorder(hx, order(setdiff(names(hx), end_cols)))

    fwrite(hx, paste0(dir, "/history.csv"), row.names = TRUE)

    cat("\n")

    # check for improvement
    better <- hx[[iteration, "score"]] == max(hx[["score"]])

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
  setcolorder(history, order(setdiff(names(history), c(end_cols, cwp))))
  best_cwp <- history[score == max(score), cwp][[1]]
  list(history = history, best_cwp = best_cwp)

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
    prefilter_int      = 10000,  # originally 100
    mzdiff             = c(-0.001, 0.01),
    noise              = 0  # originally 0, increase to increase speed
  )
}


param_types_centwave <- function() {

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

  param_types_centwave()

  # check parameter names
  nms <- names(parameter_list)
  bad_names <- nms[nms %nin% c(quant, qual, lists)]
  if (length(bad_names) != 0) {
    stop("These parameters are not used by centWave: ",
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
  for (nm in nms){
    if (nm %in% c("ppm")) {
      if (min(parameter_list[[nm]]) <= 0) {
        stop(paste("The parameter", nm, "must be greater than 0"))
      }
    } else if (nm %in% c("min_peakwidth", "max_peakwidth")) {
      if (min(parameter_list[[nm]]) <= 0) {
        stop(paste("The parameter", nm, "must be greater than 0"))
      }
    } else if (nm %in% c("snthresh", "prefilter_k", "prefilter_int", "noise")) {
      if (min(parameter_list[[nm]]) < 0) {
        stop(paste("The parameter", nm, "must not be negative"))
      }
    }
  }
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
