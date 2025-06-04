#' This function normalises each column in F to row sum
#' 
#' @keywords internal
#' 
#' @param Fmat A matrix or data.frame with the last column containing non-zeros
#'
#' @return A list consisting of two components:
#'     - a matrix of pigment ratios normalized to row sums
#'     - a vector of row sums
#'
#' @examples
Normalise_F <- function(Fmat) {
  Fmat  <- as.matrix(Fmat)           # convert to matrix
  F_1   <- Fmat / Fmat[, ncol(Fmat)] # divide Fmat by last column
  F.sum <- rowSums(F_1)              # sum rows
  F_1   <- F_1 / F.sum               # divide by sum
  return(list(F_1, F.sum))
}
