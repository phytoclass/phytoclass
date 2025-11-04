#' Apply the steepest descent algorithm to optimize pigment ratios
#' in phytoplankton classification. Loosely wraps Seepest_Descent function.
#'  
#' @keywords internal  
#'
#' @param Ft Initial F matrix containing pigment ratios
#' @param min.value UNUSED?
#' @param max.value UNUSED?
#' @param place Vector of indices where F matrix has non-zero values
#' @param S Matrix of sample measurements
#' @param cm Vector of bounded weights for each pigment
#' @param num.loops Maximum number of iterations for the steepest descent
#'
#' @return A list containing:
#'   \code{1}: The optimized F matrix
#'   \code{2}: Final RMSE value
#'
#' @examples
#'  Fmat <- as.matrix(phytoclass::Fm)
#'  S <- as.matrix(phytoclass::Sm)
#'  S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
#'  place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
#'  num.loops <- 2
#'  # Run SAALS
#'  result <- phytoclass:::SAALS(Fmat, NULL, NULL, place, S, S_weights, num.loops)

SAALS <- function(Ft, min.value, max.value, place, S, cm, num.loops){
  g <- Steepest_Descent(Ft, place, S, cm, num.loops)
  err <- g[[2]]
  g <- g[[1]]
  #if (length(d) >0){
  #  g <- Ft
  #}
  return(list(g,err))
}
