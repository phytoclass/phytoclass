#' Apply the steepest descent algorithm to optimize pigment ratios
#' in phytoplankton classification
#'  
#' @keywords internal  
#'
#' @param Ft Initial F matrix containing pigment ratios
#' @param min.value Minimum allowed values for pigment ratios
#' @param max.value Maximum allowed values for pigment ratios
#' @param place Vector of indices where F matrix has non-zero values
#' @param S Matrix of sample measurements
#' @param cm Vector of bounded weights for each pigment
#' @param num.loops Maximum number of iterations for the steepest descent
#'
#' @return A list containing:
#'   \item{[[1]]}{The optimized F matrix}
#'   \item{[[2]]}{Vector of error values over iterations}
#'
#' @examples
#' # Create sample matrices
#' Ft <- matrix(c(0.5, 0.3, 0.2,
#'                0.4, 0.1, 0.5), nrow=2, byrow=TRUE)
#' min.value <- matrix(0.1, nrow=2, ncol=3)
#' max.value <- matrix(0.9, nrow=2, ncol=3)
#' place <- c(1, 2, 3, 4, 5, 6)  # all elements
#' S <- matrix(runif(12), nrow=4)  # 4 samples, 3 pigments
#' cm <- c(1, 1, 1)  # equal weights
#' 
#' # Run optimization
#' result <- SAALS(Ft, min.value, max.value, place, S, cm, num.loops=100)

SAALS <- function(Ft, min.value, max.value, place, S, cm, num.loops){
  g <- Steepest_Descent(Ft, place, S, cm, num.loops)
  err <- g[[2]]
  g <- g[[1]]
  #if (length(d) >0){
  #  g <- Ft
  #}
  return(list(g,err))
}
