#' Conduit between minimise_elements function and Fac_F_R 
#' of steepest descent algorithm.
#' @param F     
#' @param place   
#' @param S   
#' @param cm   
#'
#' @return
#' @export
#'
#' @examples
Conduit_2 <- function(F, place, S, cm){
  F.locs <- vector()
  F.old <- NNLS_MF(F, S, cm)
  F.news <- Fac_F_RR2(F.old, vary = place, place, S, cm)
  F.new <- F.news[[1]]
  n <- F.news[[2]]
  res <- list(F.new, n, F.old)
  return(res)
}
