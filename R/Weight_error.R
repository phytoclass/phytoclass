#' Add weights
#'
#' @param S 
#'
#' @return A matrix
#' @export
#'
#' @examples
Weight_error <- function(S){
  S <- S%*%diag(cm)
  return(S)
}