#' Apply randomization functions to matrices
#'
#' @param F
#' @param i 
#' @param S
#' @param cm
#'
#' @return
#' @export
#'
#' @examples
Replace_Rand1 <- function(F, i, S, cm){
  # randomise first element of matrix  
  F.new <- as.matrix(replace(F[[1]], i, rand(F[[1]][i], min.scaler = 0.99, max.scaler = 1.01))) 
  F.new <- Fac_F(F.new, S, cm)
  # Which elements decrease the error? Store the location of the elements that decrease it   
  v <- which(F.new[[2]] < F[[2]])
  res <- c(F.new,v) 
  return(res)
}
