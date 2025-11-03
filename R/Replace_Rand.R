#' Select the new F matrix element with lowest error in the steepest
#' descent algorithm. Randomly modifies a single element and checks if the 
#' modification reduces the error.
#' 
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
#' # Create sample matrices
#' F <- matrix(c(0.5, 0.3, 0.2,
#'               0.4, 0.1, 0.5), nrow=2, byrow=TRUE)
#' S <- matrix(runif(12), nrow=4)
#' cm <- c(1, 1, 1)
#' Fmat <- list(F, 0.1, matrix(1, nrow=2, ncol=4))
#' 
#' # Try modifying element at index 1
#' result <- Replace_Rand(Fmat, 1, S, cm, 0.99, 1.01)
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

#' Randomise individual elements in the F matrix by applying scaling factors
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
#' new_value <- Randomise_elements(x, 0.99, 1.01)  # +/- 1% change
#' 
#' # Handle small values
#' small_x <- 0.0005
#' new_small <- Randomise_elements(small_x, 0.99, 1.01)  # will use 0.001
Randomise_elements <- function(x, min.scaler, max.scaler) {
  x[x < 0.001] <- 0.001
  round(runif(n = 1, min = x * min.scaler, max = x * max.scaler), digits = 4)
}
