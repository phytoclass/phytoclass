#' Title
#' 
#' @keywords internal
#'
#' @param x   xx
#'
#' @return
#'
#' @examples
target = function(x){
  return(ifelse(x<0,0,exp(-x)))
}
