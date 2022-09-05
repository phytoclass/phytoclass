#' Same as before, but without the 'N' input
#'
#' @param Fn 
#' @param Temp 
#' @param chlv 
#' @param s_c 
#' @param S  
#' @param cm
#' @param min.val 
#' @param max.val
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom stats runif
Random_neighbour2 <- function(Fn, Temp, chlv, s_c, place, S, cm, min.val, max.val){
  s_c <- vectorise(s_c[,1:ncol(s_c)-1])
  #### vectorise function outputs all non-zero elements as a vector (excluding chl column)  
  SE <-   Wrangling(Fn, min.val, max.val)[[3]] 
  minF <- Wrangling(Fn, min.val, max.val)[[1]]
  maxF <- Wrangling(Fn, min.val, max.val)[[2]]
  ki <- maxF-minF
  rand <- round(runif(n = length(s_c), -1, 1),4)
  SA <- (SE + (Temp*0.5) *ki*rand)
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
    SA2 <- (SE[d] +(Temp*0.5)*kir*nr)
    SA[d] <- SA2 
    d <- which(SA < minF | SA > maxF)
    #print(loop)
    if (loop > 2000){
      nn <- (minF[d]+maxF[d])/2
      f <- round(runif(n=length(d),minF[d],maxF[d]),4)
      SA[d] <- f
      d <- which(SA < minF | SA > maxF)
    }
  }
  Fn <- Fn[,1:ncol(Fn)-1] #### If error is lower, reassign the values
  Fn[Fn >0] <- SA
  Fn <- cbind(Fn,chlv)
  colnames(Fn) <- colnames(F)
  F.n <- Fac_F(Fn, S, cm)
  return(F.n, S, cm)
}
