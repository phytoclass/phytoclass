#' Part of the steepest descent algorithm and work to reduce error given 
#' the S and F matrices.
#' 
#' @keywords internal 
#'
#' @param Fmat XX
#' @param vary XX
#' @param S XX
#' @param cm XX
#'
#' @return
#'
#' @examples
# Inputs are F and which elements to vary, should be all elements
Fac_F_RR1 <- function(Fmat, vary, S, cm){ 
  F.locs <- vector() 
  Fs <- Fmat[[1]]
  F.new <- lapply(vary, function(i){ # randomises every element in 'vary' 
    Replace_Rand(Fmat, i, S, cm, min.scaler = 0.99, max.scaler = 1.01)})
  # Shows which elements reduce error  
  cont <- lapply(1:length(F.new), function(i){ c <- which(length(F.new[[i]]) == 4)}) 
  conts <- which(cont==1)
  # Procedure for if no elements reduce error  
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
      F.new <- NNLS_MF(Fmat[[1]], S, cm)
      cont <- vary
    }
  }
  else{
    F.new <- NNLS_MF(Fmat[[1]], S, cm)
    cont <- vary
  }
  res <- list(F.new,cont)
  return(res)
  
}

