#' Select the new F matrix element with lowest error in the steepest
#' descent algorithm. 
#' 
#' @keywords internal
#'
#' @param Fmat   xx
#' @param i xx
#' @param S   xx
#' @param cm   xx
#' @param min.scaler   xx
#' @param max.scaler  xx
#'
#' @return
#'
#' @examples
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