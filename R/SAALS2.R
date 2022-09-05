#' to apply to SD algorithm after each iteration of simulated annealing.
#'
#' @param Ft 
#'
#' @return
#' @export
#'
#' @examples
SAALS2 <- function(Ft, S, cm){
  g <- Try_This2(Ft, F, cm)
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
