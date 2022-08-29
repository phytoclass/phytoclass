#' Apply randomisation function to matrices and keep element that reduce error x4
#'
#' @param F
#' @param vary 
#' @param S
#' @param cm
#'
#' @return
#' @export
#'
#' @examples
Fac_F_RR1 <- function(F, vary, S, cm){ # Inputs are #F and which elements to vary, should be all elements
  F.locs <- vector() 
  Fs <- F[[1]]
  F.new <- lapply(vary, function(i){ # randomises every element in 'vary' 
    Replace_Rand1(F,i)})
  cont <- lapply(1:length(F.new), function(i){ c <- which(length(F.new[[i]]) == 4)}) # Shows which elements reduce error
  conts <- which(cont==1)
  if(!is.null(length(conts))){ # Procedure for if no elements reduce error
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
      F.new <- Fac_F(F.new, S, cm)
    }
    else{
      F.new <- Fac_F(F[[1]], S, cm)
      cont <- vary
    }
  }
  else{
    F.new <- Fac_F(F[[1]], S, cm)
    cont <- vary
  }
  res <- list(F.new,cont)
  return(res)
  
}

