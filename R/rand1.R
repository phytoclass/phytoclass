#' Randomisation function
#'
#' @param x 
#' @param na.rm 
#'
#' @return numeric
#' @export
#'
#' @examples
rand1 <- function(x, na.rm = FALSE){x <- ifelse(x<0.001,0.001,x)
round(runif(n=1, x*0.99,x*1.01),4)}