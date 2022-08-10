
# Condisder a random neighbour to the inital value by randomising from normal distrubtion with scaling factors etc... 
# This has the input N and is called in the simuated annealing output to randomise an element that is outside its min/max value
Random_neighbour <- function(Fn,Temp,chlv,s_c,N){
  k <- match(N,place)
  s_c <- s_c[N]
  SE <- Wrangling(Fn)[[3]] #### vectorise function outputs all non-zero elements as a vector (excluding chl column)
  SE <- SE[k]
  minF <- Wrangling(Fn)[[1]]
  minF <- minF[k]
  maxF <- Wrangling(Fn)[[2]]
  maxF <- maxF[k]
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
    #print(loop) # if loop is > 3000, just select it from a uniform distribution between the max and min
    if (loop > 3000){
      nn <- minF[d]+maxF[d]/2
      f <- round(runif(n=length(d),nn*0.01,nn*1.5),4)
      SA[d] <- f
      d <- which(SA < minF | SA > maxF)
    }
  }
  Fn <- Fn[,1:ncol(Fn)-1] #### If error is lower, reassign the values
  Fn[N] <- SA
  Fn <- cbind(Fn,chlv)
  colnames(Fn) <- colnames(F)
  F.n <- Fac_F(Fn)
  return(F.n)
}


