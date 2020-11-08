score_peaks <- function(xcmsnexp) {
  peak_df <- data.table(xcms::chromPeaks(xcmsnexp), keep.rownames = TRUE)

  a <-
    peak_df[

      # summarize peak tables by sample
      , list(peak_tables = list(.SD)), by = sample
    ][

      # count total number of peaks
      , total := purrr::map_dbl(peak_tables, nrow)
    ][

      # identify isotope peaks
      , isotope_list := purrr::map(peak_tables, find_isotopes)
    ][

      # identify low intensity peaks with likely undetectable isotopes
      , indeterminate := purrr::map2_dbl(peak_tables, isotope_list, find_indeterminate_peaks)
    ][

      # count isotopes in list
      , isotopes := purrr::map_dbl(isotope_list, length)

    ]

  a[
    , lapply(.SD, sum), .SDcols = c("total", "indeterminate", "isotopes")
  ][
    , score := isotopes ^ 2 / (total - indeterminate)
  ]
}


find_isotopes <- function(peak_table) {

  peak_table <- peak_table[order(peak_table$mz)]

  pair <- list()
  counter <- 0

  C13_mass_defect <- 1.003355
  ppm <- 5

  for (row in 1:nrow(peak_table)) {

    query <- peak_table[row]

    mz <- peak_table[[row, "mz"]]

    mz_low <- mz + C13_mass_defect - mz * ppm / 1E6
    mz_high <- mz + C13_mass_defect + mz * ppm / 1E6

    matches <- peak_table[mz >= mz_low & mz <= mz_high][, id := .I]

    if (nrow(matches) > 0) {
      counter <- counter + 1
      query_2 <- query[rep(1, nrow(matches))][, id := .I]
      pair[[counter]] <- merge(
        query_2,
        matches,
        suffixes = c("_1", "_2"),
        by = "id"
      )
    }

  }

  isotope_pairs <- data.table::rbindlist(pair)

  # identify overlapping retention time windows
  rt_overlap <-
    pmax(isotope_pairs$rtmin_1, isotope_pairs$rtmin_2) <
    pmin(isotope_pairs$rtmax_1, isotope_pairs$rtmax_2)

  # confirm expected peak intensity ratios
  intb_ratio <- isotope_pairs$intb_2 / isotope_pairs$intb_1
  carbon <- maximize_carbons(isotope_pairs$mz_1)
  predicted_ratio <- predict_isotope_ratio(carbon)
  expected_ratio <-
    intb_ratio >= 0.75 * predict_isotope_ratio(1) &
    intb_ratio <= 1.25 * predicted_ratio

  # list unique peak ids
  unique(
    unlist(
      isotope_pairs[rt_overlap & expected_ratio, c("rn_1", "rn_2")]
    ),
    use.names = FALSE
  )
}


maximize_carbons <- function(mass) {
  C <- 12.0
  H  <- 1.0078250170
  CH3 <- C + H * 3
  CH2 <- C + H * 2
  floor((mass - 2 * CH3) / CH2) + 2
}


predict_isotope_ratio <- function(carbons) {
  C13_abundance <- 0.0107
  dbinom(1, carbons, C13_abundance) / dbinom(0, carbons, C13_abundance)
}


find_indeterminate_peaks <- function(peak_tbl, isotope_list) {

  no_isotopes <- peak_tbl[peak_tbl$rn %nin% isotope_list, ]

  three_percent <- round(nrow(no_isotopes) * 0.03)

  threshold <- no_isotopes[order(intb)][1:three_percent, mean(intb)]

  carbon <- maximize_carbons(no_isotopes$mz)
  predicted_ratio <- predict_isotope_ratio(carbon)
  predicted_intb <- no_isotopes$intb * predicted_ratio
  nrow(no_isotopes[predicted_intb < threshold, ])

}
