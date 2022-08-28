#' conduit between minimise_elements function and Fac_F_R
#'
#' @return
#' @export
#'
#' @examples
Test3 <- function(F, place){
  F.old <- Fac_F(F)
  F.news <- Fac_F_RR3(F.old, place)
  F.new <- F.news[[1]]
  n <- F.news[[2]]
  res <- list(F.new, n, F.old)
  return(res)
}
