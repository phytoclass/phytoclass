#' add weights to the data.
#'
#' @param S 
#'
#' @return A matrix
#' @export
#'
#' @examples
cms<-function(S){
  n <- colMeans(S)
  S <- n^-1
  S <- ifelse(S>50,50,S)
  S[length(S)] =1
  return(S)
}