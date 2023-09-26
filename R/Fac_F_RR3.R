#' Part of the steepest descent algorithm and work to reduce error given 
#' the S and F matrices
#' 
#' @keywords internal
#'
#' @param Fmat  xx
#' @param vary xx
#' @param place   xx   
#' @param S  xx
#' @param cm xx
#'
#' @return
#'
#' @examples
Fac_F_RR3 <- function(Fmat, vary, place, S, cm){
  F.locs <- vector()
  F.new <- lapply(vary, function(i) {Replace_Rand(Fmat, i, S, cm, min.scaler = 0.97, max.scaler = 1.03)})
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
      F.new <- replace(Fmat[[1]],cont,F.news)
      F.new <- NNLS_MF(F.new, S, cm)
    }
    else{
      C <- Fac_F_RR2(Fmat, vary, place, S, cm)
      F.new <- C[[1]] 
      cont <- C[[2]]
    }
  }
  else{
    C <- Fac_F_RR1(Fmat, place)
    F.new <- C[[1]] 
    cont <- C[[2]]
  }
  res <- list(F.new,cont)
  return(res)
  
}
