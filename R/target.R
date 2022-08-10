#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
target = function(x){
  return(ifelse(x<0,0,exp(-x)))
}
