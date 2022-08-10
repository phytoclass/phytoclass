#' Apply randomisation function to matrices and keep element that reduce error x4
#'
#' @param vary 
#'
#' @return
#' @export
#'
#' @examples
Fac_F_RR4 <- function(F,vary){
  F.locs <- vector()
  Fs <- F[[1]]
  #  if (is.null(vary)){vary <- place}
  F.new <- lapply(vary, function(i){
    Replace_Rand4(F,i)})
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
      F.new <- Fac_F(F[[1]])
      cont <- vary
    }
  }
  else{
    F.new <- Fac_F(F[[1]])
    cont <- vary
  }
  res <- list(F.new,cont)
  return(res)
  
}


