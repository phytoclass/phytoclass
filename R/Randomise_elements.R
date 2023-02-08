#' Randomisation function
#'
#' @param x 
#' @param na.rm 
#' @param min.scaler
#' @param max.scaler
#'
#' @return numeric
#' @export
#'
#' @examples
rand <- function(x, min.scaler, max.scaler){
  x <- ifelse(x<0.001, 0.001, x)
  round(runif(n=1, min = x*min.scaler, max = x*max.scaler), digits = 4)
}
