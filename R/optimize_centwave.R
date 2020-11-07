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

  # set up iteration loop
  history <- list()
  iteration <- 1
  while(iteration <= 1) {

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
        cwp, function(x) {
          xcms::findChromPeaks(
            raw_data,
            param = x,
            BPPARAM = BiocParallel::SerialParam())
        }
      )

    ## score experiment

    ## calculate model

    ## find best parameters

    ## score best parameters

    ## assign to history

    ## make sure is improved

    ## adjust intervals

    # increase counter
    iteration <- iteration + 1

  }

  # output results
  return(
    xcmsnexp
  )

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
    if (nm %in% c("ppm", "min_peakwidth", "max_peakwidth")) {
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
