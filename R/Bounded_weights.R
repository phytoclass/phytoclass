#' Add weights to the data, bound at a maximum.
#' 
#' @param S   Sample data matrix – a matrix of pigment samples
#' @param weight.upper.bound  Upper bound for weights (default is 30)        
#'
#' @return A vector with upper bounds for weights
#' @export
#'
#' @examples
#' Bounded_weights(Sm, weight.upper.bound = 30)
#' 
Bounded_weights <-function(S, weight.upper.bound = 30){
  n <- colMeans(S)
  S <- n^-1
  S <- ifelse(S > weight.upper.bound, weight.upper.bound, S)
  S[length(S)] = 1
  return(S)
}
