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


#' Prochloro Wrangling
#' 
#' @keywords internal
#'
#' @param Fl 
#' @param min.val 
#' @param max.val 
#'
#' @return
#'
#' @examples
Prochloro_Wrangling <- function(Fl, min.val, max.val) {
  
  # set up initial F, Fmin, Fmax, matrix by removing Tchla column
  Fd <- Fl
  Fmin <- Fmax <- as.matrix(Fl)[, -ncol(Fl)]
  
  Fmin[Fmin > 0] <- min.val
  Fmax[Fmax > 0] <- max.val
  
  chlv  <- Fl[, ncol(Fl)]      # Chl a
  chlvp <- Fl[, ncol(Fl) - 1]  # dvChl a
  chlep <- Fl[, ncol(Fl) - 1]
  chlep[length(chlep)] <- 1
  
  # Non-Prochloro rows scaled by Chl a (exclude dvChl column)
  Fmin[-nrow(Fmin), -ncol(Fmin)] <- Fmin[-nrow(Fmin), -ncol(Fmin)] * chlv[-length(chlv)]
  Fmax[-nrow(Fmax), -ncol(Fmax)] <- Fmax[-nrow(Fmax), -ncol(Fmax)] * chlv[-length(chlv)]
  
  # Prochlorococcus row scaled by dvChl a
  Fmin[nrow(Fmin), ] <- Fmin[nrow(Fmin), ] * chlvp[length(chlvp)]
  Fmax[nrow(Fmax), ] <- Fmax[nrow(Fmax), ] * chlvp[length(chlvp)]
  
  # Vectorised CORE (base-R column-major) and remove zeros
  Fmin_vec <- vectorise(Fmin)
  Fmax_vec <- vectorise(Fmax)
  SE_vec   <- vectorise(Fl[, 1:(ncol(Fmin) - 2)])
  
  return(list(Fmin_vec, Fmax_vec, SE_vec, chlv, chlvp))
}