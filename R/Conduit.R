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
Conduit <- function(Fmat, place, S, cm, c_num = c(1, 2, 3)) {
  F.old <- NNLS_MF(Fmat, S, cm)
  place1 <- NULL
  if (c_num != 1) {place1 <- place }
  
  F.news <- Fac_F_RR(F.old, vary = place, place = place1, S, cm, fac_rr = c_num)
  
  res <- list(F.news[[1]], F.news[[2]], F.old)
  return(res)
}