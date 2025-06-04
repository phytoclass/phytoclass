#' Normalise matrix to row sum
#' 
#' This function normalises each column in S to row sum
#' 
#'
#' @keywords internal
#'
#' @param S   A matrix or data.frame
#'
#' @return A matrix
#'
#' @examples
#' 
Normalise_S <- function(S){
  # Normalise to unit row sum
  S <- as.matrix(S)
  S <- S / rowSums(S)
  return(S)
}


# Normalise_S2 <- function(S){
#   S <- S[,1:ncol(S)]
#   #Normalise to unit row sum
#   S.sum <- rowSums(S)
#   S <- S/S.sum
#   #vectorise / matrices
#   Sm <- as.matrix(S)
#   Sn <- as.vector(Sm[1,])
#   return(Sm)
# }

