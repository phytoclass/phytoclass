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
Replace_Rand3 <- function(F, i, S, cm){
  # randomise first element of matrix  
  F.new <- as.matrix(replace(F[[1]], i, rand(F[[1]][i], min.scaler = 0.97, max.scaler = 1.03))) 
  F.new <- Fac_F(F.new, S, cm)
  v <- which(F.new[[2]] < F[[2]])
  res <- c(F.new,v)
  return(res) 
}
