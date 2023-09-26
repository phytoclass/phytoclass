#' Conduit between minimise_elements function and Fac_F_R 
#' of steepest descent algorithm.
#' 
#' @keywords internal
#'
#' @param Fmat  xx   
#' @param place  xx
#' @param S  xx
#' @param cm  xx
#' @return
#'
#' @examples
Conduit_1 <- function(Fmat, place, S, cm){
  F.locs <- vector()
  F.old <- NNLS_MF(Fmat, S, cm)
  F.news <- Fac_F_RR1(F.old, place, S, cm)
  F.new <- F.news[[1]]
  n <- F.news[[2]]
  res <- list(F.new,n, F.old)
  return(res)
}

