score_align_group <- function(xcmsnexp) {
  f_defs <- xcms::featureDefinitions(xcmsnexp)

  # Initialize all scores
  rcs <- 0
  good <- 0
  bad <- 0
  gs <- 0

  if (nrow(f_defs) > 0) {
    peaks <- xcms::chromPeaks(xcmsnexp)

    # score rt deviation
    rt_shift <-
      lapply(f_defs$peakidx, function(x) {
        sum(abs(stats::median(peaks[x, "rt"]) - peaks[x, "rt"])) / length(x)
      })

    if (sum(unlist(rt_shift)) != 0) {
      rcs <- length(rt_shift) / sum(unlist(rt_shift))
    } else {
      rcs <- 0  # avoid division by zero
    }

    # score grouping quality
    f_summary <- xcms::featureSummary(xcmsnexp)
    good <- sum(f_summary[, "perc"] == 100 & f_summary[, "multi_perc"] == 0)
    bad <- nrow(f_summary) - good
    gs <- if (bad != 0) good^2 / bad else 0
  }

  data.table::data.table(rcs = rcs, good = good, bad = bad, gs = gs)
}
