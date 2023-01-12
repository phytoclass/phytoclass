#' Apply randomization functions to matrices
#'
#' @param F
#' @param i 
#' @param S
#' @param cm
#' @param min.scaler
#' @param max.scaler
#'
#' @return
#' @export
#'
#' @examples
Replace_Rand <- function(F, i, S, cm, min.scaler, max.scaler){
  # randomise first element of matrix  
  F.new <- as.matrix(replace(F[[1]], i, rand(F[[1]][i], min.scaler, max.scaler))) 
  F.new <- Fac_F(F.new, S, cm)
  # Which elements decrease the error? Store the location of the elements that decrease it   
  v <- which(F.new[[2]] < F[[2]])
  res <- c(F.new,v) 
  return(res)
}
