#' Select the new F matrix element with lowest error in the steepest
#' descent algorithm. 
#'
#' @param F   xx
#' @param i xx
#' @param S   xx
#' @param cm   xx
#' @param min.scaler   xx
#' @param max.scaler  xx
#'
#' @return
#'
#' @examples
Replace_Rand <- function(F, i, S, cm, min.scaler, max.scaler){
  # randomise first element of matrix  
  F.new <- as.matrix(replace(F[[1]], i, Randomise_elements(F[[1]][i], min.scaler, max.scaler))) 
  F.new <- NNLS_MF(F.new, S, cm)
  # Which elements decrease the error? Store the location of the elements that decrease it   
  v <- which(F.new[[2]] < F[[2]])
  res <- c(F.new,v) 
  return(res)
}
