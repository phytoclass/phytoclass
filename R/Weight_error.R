#' Apply weights to F/S matrices
#'
#' @param S  xx
#' @param cm  xx
#'
#' @return A matrix
#' @export
#'
#' @examples
Weight_error <- function(S, cm){
  S <- S%*%diag(cm)
  return(S)
}