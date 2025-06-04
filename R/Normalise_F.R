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


# Normalise_F2 <- function(Fmat){
#   F_1<- Fmat
#   Fchl <- F_1[,ncol(F_1)]
#   F_1 <- F_1 / Fchl
#   F.sum <- rowSums(F_1)
#   #F_1 <- F/norm(F,'F')
#   F_1 <- F_1 / F.sum
#   F_1m <- as.matrix(F_1)
#   return(list(F_1m,F.sum))
# }
