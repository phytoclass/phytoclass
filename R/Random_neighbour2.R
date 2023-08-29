#' Selects a random neighbour for the simulated annealing algorithm.
#' 
#' @keywords internal
#'
#' @param Fn   xx
#' @param Temp xx
#' @param chlv xx
#' @param s_c xx
#' @param place  xx
#' @param S      xx
#' @param cm   xx
#' @param min.val    xx
#' @param max.val   xx
#'
#' @return
#'
#' @examples
#' @importFrom stats runif
Random_neighbour2 <- function(Fn, Temp, chlv, s_c, place, S, cm, min.val, max.val){
  s_c <- vectorise(s_c[,1:ncol(s_c)-1])
  SE <-   Wrangling(Fn, min.val, max.val)[[3]] 
  minF <- Wrangling(Fn, min.val, max.val)[[1]]
  maxF <- Wrangling(Fn, min.val, max.val)[[2]]
  ki <- maxF-minF
  rand <- round(runif(n = length(s_c), -1, 1),4)
  SA <- (SE + (Temp) *ki*rand)
  SA <- as.vector(unlist(SA))
  d <- which(SA < minF | SA > maxF)
  length(d)
  loop <- 1
  while (length(d) >0){
    loop = loop +1
    nr <- round(runif(length(d),-1,1),4)
    minr <- as.vector(unlist(minF))
    maxr <- as.vector(unlist(maxF))
    mind <- minr[d]
    maxd <- maxr[d]
    kir <- maxd-mind
    SA2 <- (SE[d] +(Temp)*kir*nr)
    SA[d] <- SA2 
    d <- which(SA < minF | SA > maxF)
    #print(loop)
    if (loop > 50){
      nn <- (minF[d]+maxF[d])/2
      f <- round(runif(n=length(d),(minF[d]*1.20),(maxF[d]*0.80)),4)
      SA[d] <- f
      d <- which(SA < minF | SA > maxF)
    }
  }
  Fn <- Fn[,1:ncol(Fn)-1] #### If error is lower, reassign the values
  Fn[Fn >0] <- SA
  Fn <- cbind(Fn,chlv)
  colnames(Fn) <- colnames(F)
  F.n <- NNLS_MF(Fn, S, cm)
  return(F.n)
}

