#' Vectorise a matrix and keep non-zero elements
#' 
#' Turn each non-zero element of the F-matrix into a vector
#' 
#' @keywords internal
#'
#' @param Fmat  A matrix to vectorise
#'
#' @return A vector of non-zero pigment elements
#'
#' @examples
vectorise <- function(Fmat) {
  return(Fmat[Fmat > 0])
}
