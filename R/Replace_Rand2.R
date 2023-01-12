#' Apply randomization functions to matrices
#'
#' @param i 
#'
#' @return
#' @export
#'
#' @examples
Replace_Rand2 <- function(F, i, S, cm){
  # randomise first element of matrix  
  F.new <- as.matrix(replace(F[[1]], i, rand(F[[1]][i], min.scaler = 0.98, max.scaler = 1.02))) 
  F.new <- Fac_F(F.new, S, cm)
  v <- which(F.new[[2]] < F[[2]])
  res <- c(F.new,v)
  return(res)
}