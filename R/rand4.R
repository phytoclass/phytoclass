#' Randomisation function
#'
#' @param x 
#' @param na.rm 
#'
#' @return numeric
#' @export
#'
#' @examples
rand4 <-  function(x, na.rm = FALSE){x <- ifelse(x<0.001,0.001,x)
round(runif(n=1, x*.9,x*1.1),4)}