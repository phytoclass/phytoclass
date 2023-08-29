#' Apply the steepest descent algorithm 
#'  
#' @keywords internal  
#'
#' @param Ft xx
#' @param min.value   xx
#' @param max.value   xx
#' @param place   xx
#' @param S  xx
#' @param cm   xx
#' @param num.loops xx
#'
#' @return
#'
#' @examples


SAALS <- function(Ft, min.value, max.value, place, S, cm, num.loops){
  g <- Steepest_Descent(Ft, place, S, cm, num.loops)
  err <- g[[2]]
  g <- g[[1]]
  #if (length(d) >0){
  #  g <- Ft
  #}
  return(list(g,err))
}
