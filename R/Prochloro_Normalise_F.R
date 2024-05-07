#' Normalise F for prochloro
#'
#' @keywords internal
#'
#' @param Fmat 
#'
#' @return
#'
#' @examples
Prochloro_Normalise_F <- function (Fmat) 
{
  F_1 <- Fmat
  Fchl <- F_1[, ncol(F_1)]
  F_1[1:nrow(F_1)-1,] <- F_1[1:nrow(F_1)-1,]/Fchl[length(Fchl)-1]
  F.sum <- rowSums(F_1)
  F_1 <- F_1/F.sum
  F_1m <- as.matrix(F_1)
  return(list(F_1m, F.sum))
}

