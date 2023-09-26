#' This function normalises each column in F to row sum
#' 
#' @keywords internal
#' 
#' @param Fmat     xx
#'
#' @return A matrix
#'
#' @examples
Normalise_F <- function(Fmat){
  F_1<- Fmat
  Fchl <- F_1[,ncol(F_1)]
  F_1 <- F_1 / Fchl
  F.sum <- rowSums(F_1)
  #F_1 <- F/norm(F,'F')
  F_1 <- F_1 / F.sum
  F_1m <- as.matrix(F_1)
  return(list(F_1m,F.sum))
}
