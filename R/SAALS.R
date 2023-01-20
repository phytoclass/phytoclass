#' to apply to SD algorithm after each iteration of simulated annealing.
#'

#' @param Ft 
#' @param min.value
#' @param max.value
#' @param place
#' @param S  
#' @param cm
#' @param num.loops
#'
#' @return
#' @export
#'
#' @examples


SAALS <- function(Ft, min.value, max.value, place, S, cm, num.loops){
  g <- Steepest_Descent (Ft, place, S, cm, num.loops)
  err <- g[[2]]
  g <- g[[1]]
  gchl <- g[,ncol(g)]
  gn <- g / gchl
  n <- vectorise(gn)
  d <- which(n < min.value | n > max.value)
  #if (length(d) >0){
  #  g <- Ft
  #}
  return(list(g,err))
}
