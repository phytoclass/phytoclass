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
  F.new <- as.matrix(replace(F[[1]],i,rand1(F[[1]][i]))) # randomise first element of matrix
  F.new <- Fac_F(F.new, S, cm)
  v <- which(F.new[[2]] < F[[2]]) #Which elements decrease the error? Store the location of the elements that decrease it 
  res <- c(F.new,v) 
  return(res)
}
