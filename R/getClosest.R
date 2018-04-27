#' Finds the Closest Genetic Match for Each Sample
#'
getClosest <- function(i) {
  which(i == max(i, na.rm = TRUE))
}
