#' A function that reduces every for every element that didn't reduce in index function
#'
#' @keywords internal
#'
#' @param Fmat   xx
#' @param place  xx
#' @param S    xx
#' @param cm     xx 
#'
#' @return
#'
#' @examples
Minimise_elements1 <- function(Fmat, place, S, cm){  
  # A function that reduces every for every element that didn't reduce in index function
  f <- Conduit_1(Fmat, place, S, cm) # Calls index function
  F.new <- f[[1]] # F matrix
  n <- f[[2]] #elements that reduce error
  if (is.null(n)){n <- place}
  F.old <- f[[3]] # old F matrix
  F.initial <- F.new # Fac_F new matrix
  # Fac_F new matrix
  g <-Fac_F_RR1(F.new, place, S, cm)
  if (g[[1]][[2]] < F.initial[[2]]){F.new <- g[[1]]}
  n <- g[[2]]
  res <- list(F.new,n)
  return(F.new)
}
