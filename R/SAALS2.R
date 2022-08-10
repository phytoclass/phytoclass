
SAALS2 <- function(Ft){
  g <- Try_This2(Ft)
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
