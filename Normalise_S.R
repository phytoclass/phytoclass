
#' This function normalises each column in S to row sum
#' 
#'
#' @param S 
#'
#' @return A matrix
#' @export
#'
#' @examples
#' 
Normalise_S <- function(S){
  S <- S[,1:ncol(S)]
  #Normalise to unit row sum
  S.sum <- rowSums(S)
  S <- S/S.sum  
  #vectorise / matrices
  Sm <- as.matrix(S)
  Sn <- as.vector(Sm[1,])
  return(Sm)
}
