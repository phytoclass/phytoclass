#' conduit between minimise_elements function and Fac_F_R
#'
#' @return
#' @export
#'
#' @examples
Test4 <- function(F, NSA){
  F.old <- Fac_F(F)
  F.news <- Fac_F_RR4(F.old, NSA)
  F.new <- F.news[[1]]
  n <- F.news[[2]]
  res <- list(F.new, NSA, F.old)
  return(res)
}
