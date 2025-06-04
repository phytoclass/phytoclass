#' Wrangle data to vectors
#' 
#' Converts data-types and selects data for randomisation in 
#' the simulated annealing algorithm
#' 
#' @keywords internal
#'
#' @param Fl      A matrix of the initial F matrix (i.e. pigment ratio matrix)
#' @param min.val A vector of the minimum values for each non-zero pigment ratios 
#' @param max.val A vector of the maximum values for each non-zero pigment ratios
#'
#' @return
#'     A list containing following components:
#'     - A vector Fmin with the minimum pigment ratio values
#'     - A vector Fmax with the maximum pigment ratio values
#'     - A vector SE with the current pigment ratio values
#'     - A vector chlv with the pigment ratio values for the last column in Fl
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
