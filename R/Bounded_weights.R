#' Add weights to the data, bound at a maximum.
#'
#' @param S XX         
#' @param weight.upper.bound  XX        
#'
#' @return A matrix
#'
#' @examples
Bounded_weights <-function(S, weight.upper.bound){
  n <- colMeans(S)
  S <- n^-1
  S <- ifelse(S > weight.upper.bound, weight.upper.bound, S)
  S[length(S)] = 1
  return(S)
}
