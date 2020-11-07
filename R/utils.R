#' Negate %in%
#'
#' This function will return a logical vector that is the negation of %in%.
#'
#' @param x A vector
#' @param table A vector matching the mode of `x`
#'
#' @return A logical vector with length = `length(x)`
#' @keywords internal
#'
#' @examples
#' c(1, 2, 3) %nin% c(3, 4, 5)
#'
"%nin%" <- function(x, table) {
  match(x, table, nomatch = 0L) == 0L
}
