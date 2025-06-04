#' Converts data-types and selects data for randomisation in 
#' the simulated annealing algorithm
#' 
#' @keywords internal
#'
#' @param Fl      The initial F matrix (i.e. pigment ratio matrix)
#' @param min.val The minimum values for each pigment ratio 
#' @param max.val The maximum values for each pigment ratio 
#'
#' @return
#'
#' @examples
Wrangling <- function(Fl, min.val, max.val) {
  
  # set up initial F, Fmin, Fmax, matrix by removing Tchla column
  Fd <- Fmin <- Fmax <- as.matrix(Fl)[, -ncol(Fl)]
  
  # set all non-zero elements of F matrix to the minimum and maximum values
  Fmin[Fmin > 0] <- min.val
  Fmax[Fmax > 0] <- max.val
  
  # extract chlorophyll-a values once weighted to rowsums for initial F matrix
  chlv <- Fl[, ncol(Fl)]
  
  # multiply the minimum value by weighted chlorophyll to updated ratios
  Fmin <- Fmin * chlv
  Fmax <- Fmax * chlv
  
  # vectorise function outputs all non-zero elements as a vector (excluding chl column)
  Fmin <- vectorise(Fmin)
  Fmax <- vectorise(Fmax)
  SE   <- vectorise(Fd)
  
  return(list(Fmin, Fmax, SE, chlv))
}
