score_align_group <- function(xcmsnexp) {

  f_defs <- xcms::featureDefinitions(xcmsnexp)

  if (nrow(f_defs) == 0) {
    return(list(rcs = 0, good = 0, bad = 0, gs = 0))
  }

  peaks <- xcms::chromPeaks(xcmsnexp)

  # score rt deviation
  rt_shift <-
    lapply(f_defs$peakidx, function(x) {
      sum(abs(stats::median(peaks[x, "rt"]) - peaks[x, "rt"])) / length(x)
    })

  rcs <- length(rt_shift) / sum(unlist(rt_shift))

  # score grouping quality
  f_summary <- xcms::featureSummary(xcmsnexp)

  good <- sum(f_summary[, "perc"] == 100 & f_summary[, "multi_perc"] == 100)
  bad <- nrow(f_summary) - good

  gs <- good ^ 2 / bad

  list(rcs = rcs, good = good, bad = bad, gs = gs)

}
