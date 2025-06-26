#' Part of the steepest descent algorithm
#' 
#' @keywords internal
#'
#' @param Fmat    xx
#' @param place  xx
#' @param S   xx
#' @param cm  xx
#'
#' @return
#'
#' @examples
Minimise_elements_comb <- function(Fmat, place, S, cm, c1_num = c(1, 2, 3)) { # A function that reduces every for every element that didn't reduce in index function
  
  f     <- Conduit(Fmat, place, S, cm, c_num = c1_num) # Calls index function
  F.new <- f[[1]] # F matrix
  n     <- f[[2]] # elements that reduce error
  if (is.null(n)) {
    n <- place
  }
  
  F.old     <- f[[3]] # old F matrix
  F.initial <- F.new # Fac_F new matrix
  
  # Fac_F new matrix
  place1 <- NULL
  if (c1_num != 1) {place1 <- place }
  
  g <- Fac_F_RR(F.new, vary = place, place = place1, S, cm, fac_rr = c1_num)
  
  if (g[[1]][[2]] < F.initial[[2]]) {
    F.new <- g[[1]]
  }
  
  n   <- g[[2]]
  res <- list(F.new, n)
  return(F.new)
  
}