

# to apply to SD algorithm after each iteration of simulated annealing. 
SAALS <- function(Ft){
  g <- Try_This(Ft)
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
