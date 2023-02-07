#' Add weights to the data, bound at a maximum.
#'
#' @param S 
#' @param weight.upper.bound
#'
#' @return A matrix
#' @export
#'
#' @examples
Bounded_weights <-function(S, weight.upper.bound = 30){
  n <- colMeans(S)
  S <- n^-1
  S <- ifelse(S > weight.upper.bound, weight.upper.bound, S)
  S[length(S)] = 1
  return(S)
}
