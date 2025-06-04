#' Normalise matrix to row sum
#' 
#' This function normalises each column in S to row sum
#' 
#'
#' @keywords internal
#'
#' @param S   A matrix or data.frame
#'
#' @return A matrix
#'
#' @examples
#' 
Normalise_S <- function(S){
  # Normalise to unit row sum
  S <- as.matrix(S)
  S <- S / rowSums(S)
  return(S)
}
