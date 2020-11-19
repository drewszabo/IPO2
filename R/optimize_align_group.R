

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

    # run alignment

    # run grouping

    # score

    # calculate model

    # find maximum

    # make plots

    # assign to history

    # assess improvement

    # adjust intervals

    # increment counter
    iteration <- iteration + 1

  }

  # output results

}

get_obiwarp_defaults <- function() {
  list(
    binsize            = 1,
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
    binSize            = 0.25, # size of the overlapping slices in m/z dimension
    maxFeatures        = 50 # maximum number of peak groups to be identified in a single slice
  )
}

suggest_align_group_params <- function() {
  list(
    binsize            = c(0.7, 1),
    response           = c(1, 10),
    gapInit            = c(0, 0.4),
    gapExtend          = c(2.1, 2.7),
    bw                 = c(22, 38),
    minFraction        = c(0.3, 0.7),
    binSize            = c(0.015, 0.035)
  )
}

param_types_align_group <- function() {

  assign(
    "quant",
    c(
      "binsize",
      "response",
      "gapInit",
      "gapExtend",
      "factorDiag",
      "factorGap",
      "initPenalty",
      "bw",
      "minFraction",
      "minSamples",
      "binSize",
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
