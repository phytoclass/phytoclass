
Fac_F_RR2 <- function(F,vary){
  F.locs <- vector()
  F.new <- lapply(vary, function(i){ 
    Replace_Rand2(F,i)})
  cont <- lapply(1:length(F.new), function(i){ c <- which(length(F.new[[i]]) == 4)})
  conts <- which(cont==1)
  if(!is.null(length(conts))){
    cont <- as.list(vary[conts])
    conts <- as.list(conts)
    Locs <-cont
    F.news <- sapply(conts, function(i) {
      sapply(cont, function(j){
        F.locs[[length(F.locs)+1]] <- F.new[[i]][[1]][[j]]
      })
    })
    if (length(F.news) > 0){
      if (length(F.news) >1 ){F.news <- diag(F.news)}
      else{F.news <- F.news[[1]]} 
      cont <- unlist(cont)
      F.new <- replace(F[[1]],cont,F.news)
      F.new <- Fac_F(F.new)
    }
    else{
      C <- Fac_F_RR1(F,place)
      F.new <- C[[1]] 
      cont <- C[[2]]
    }
  }
  else{
    C <- Fac_F_RR1(F,place)
    F.new <- C[[1]] 
    cont <- C[[2]]
  }
  res <- list(F.new,cont)
  return(res)
  
}

