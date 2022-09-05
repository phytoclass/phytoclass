#' to apply to SD algorithm after each iteration of simulated annealing.
#'
#' @param Ft   
#' @param S  
#' @param Cm
#'
#' @return
#' @export
#'
#' @examples
SAALS <- function(Ft, S, cm){
  g <- Try_This(Ft, S, cm)
  err <- g[[2]]
  g <- g[[1]]
  gchl <- g[,ncol(g)]
  gn <- g / gchl
  n <- vectorise(gn)
  d <- which(n < min | n > max)
  #if (length(d) >0){
  #  g <- Ft
  #}
  return(list(g,err))
}
