#' to apply to SD algorithm after each iteration of simulated annealing.
#'

#' @param Ft 
#' @param place
#' @param min.value
#' @param max.value
#' @param S  
#' @param Cm
#'
#' @return
#' @export
#'
#' @examples


SAALS <- function(Ft, min.value, max.value, place, S, cm){
  g <- Try_This(Ft, place, S, cm, num.loops = 10)
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
