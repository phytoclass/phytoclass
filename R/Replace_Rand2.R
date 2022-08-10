#' Apply randomization functions to matrices
#'
#' @param i 
#'
#' @return
#' @export
#'
#' @examples
Replace_Rand2 <- function(F,i){
  F.new <- as.matrix(replace(F[[1]],i,rand2(F[[1]][i]))) # randomise first element of matrix
  F.new <- Fac_F(F.new)
  v <- which(F.new[[2]] < F[[2]])
  res <- c(F.new,v)
  return(res)
}