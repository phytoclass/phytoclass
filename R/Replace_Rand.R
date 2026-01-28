#' Select the new F matrix element with lowest error in the steepest
#' descent algorithm. Randomly modifies a single element and checks if the 
#' modification reduces the error.
#' 
#' @importFrom stats runif
#' @keywords internal
#'
#' @param Fmat A list containing the F matrix, RMSE, and other components
#' @param i The index of the element to modify in the F matrix
#' @param S Sample data matrix - matrix of pigment samples
#' @param cm A vector of bounded weights for each pigment
#' @param min.scaler Minimum scaling factor to apply (e.g., 0.99 for 1% decrease)
#' @param max.scaler Maximum scaling factor to apply (e.g., 1.01 for 1% increase)
#'
#' @return A list containing:
#'   \item{F matrix}{The modified F matrix}
#'   \item{RMSE}{Root mean square error of the new solution}
#'   \item{C matrix}{The concentration matrix}
#'   \item{Improved}{Logical indicating if the modification reduced error}
#'
#' @examples
#'  # Setup based on Fac_F_RR usage
#'  Fmat <- as.matrix(phytoclass::Fm)
#'  S <- as.matrix(phytoclass::Sm)
#'  cm <- as.numeric(phytoclass:::Bounded_weights(S))
#'
#'  # Get Fmat as a list from NNLS_MF (as used in Fac_F_RR)
#'  Fmat_list <- phytoclass::NNLS_MF(Fmat, S, cm)
#'
#'  # Test with a single index
#'  i <- 1 # first non-zero element to modify
#'  min.scaler <- 0.99
#'  max.scaler <- 1.01
#'
#'  # Run Replace_Rand
#'  result <- phytoclass:::Replace_Rand(Fmat_list, i, S, cm, min.scaler, max.scaler)
Replace_Rand <- function(Fmat, i, S, cm, min.scaler, max.scaler) {
  
  # randomise non-zero elements
  new_rand <- Randomise_elements(Fmat[[1]][i], min.scaler, max.scaler) # randomize one element
  F_new    <- replace(Fmat[[1]], i, new_rand)
  F_new    <- as.matrix(F_new)
  F_new    <- NNLS_MF(F_new, S, cm)
  
  # Which elements decrease the error? Store the location of the elements that decrease it
  v   <- F_new[[2]] < Fmat[[2]] # compare RMSE from previous, if better TRUE
  res <- c(F_new, v)
  
  return(res)
}

#' Randomise value by applying scaling factors
#' within specified bounds. Small values (< 0.001) are set to 0.001 to avoid
#' numerical issues.
#' 
#' @keywords internal
#'
#' @param x The element value to randomize
#' @param min.scaler Minimum scaling factor to apply (e.g., 0.99 for 1% decrease)
#' @param max.scaler Maximum scaling factor to apply (e.g., 1.01 for 1% increase)
#'
#' @return A numeric value between x*min.scaler and x*max.scaler, rounded to 4 decimals
#'
#' @examples
#' # Randomize a single value
#' x <- 0.5
#' new_value <- phytoclass:::Randomise_elements(x, 0.99, 1.01)  # +/- 1% change
#' 
#' # Handle small values
#' small_x <- 0.0005
#' new_small <- phytoclass:::Randomise_elements(small_x, 0.99, 1.01)  # will use 0.001
Randomise_elements <- function(x, min.scaler, max.scaler) {
  x[x < 0.001] <- 0.001
  round(runif(n = 1, min = x * min.scaler, max = x * max.scaler), digits = 4)
}
