#' Add weights
#'
#' @param S 
#' @param cm
#'
#' @return A matrix
#' @export
#'
#' @examples
Weight_error <- function(S, cm){
  S <- S%*%diag(cm)
  return(S)
}