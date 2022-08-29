#' Apply randomization functions to matrices
#'
#' @param i 
#'
#' @return
#' @export
#'
#' @examples
Replace_Rand3 <- function(F, i, S, cm){
  F.new <- as.matrix(replace(F[[1]],i,rand3(F[[1]][i]))) # randomise first element of matrix
  F.new <- Fac_F(F.new, S, cm)
  v <- which(F.new[[2]] < F[[2]])
  res <- c(F.new,v)
  return(res) 
}
