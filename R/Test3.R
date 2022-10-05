#' conduit between minimise_elements function and Fac_F_R
#'
#' @param F   
#' @param place 
#' @param S  
#' @param cm
#'
#' @return
#' @export
#'
#' @examples
Test3 <- function(F, place, S, cm){
  F.old <- Fac_F(F, S, cm)
  F.news <- Fac_F_RR3(F.old, place, S, cm)
  F.new <- F.news[[1]]
  n <- F.news[[2]]
  res <- list(F.new, n, F.old)
  return(res)
}
