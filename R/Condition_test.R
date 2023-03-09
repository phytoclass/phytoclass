#' Calculate the condition number ...
#'
#' @param S XX     
#' @param Fn  XX   
#' @param min.val XX     
#' @param max.val XX       
#'
#' @return
#' @export
#'
#' @examples


Condition_test <- function(S, Fn, min.val, max.val){
  
  condition_number <- function(S, Fd, min.val, max.val){
    Fz <- vectorise(Fd)
    rand <- vector()
    for (i in 1:length(Fz)){
      rand[[length(rand)+1]] <- stats::runif(1, min=min.val[i], max=max.val[i])
    }
    Fn <- Fd
    Fn[Fn >0] <- rand
    return(kappa(Fn %*% t(S)))
  }
  
  sn <- replicate(n = 10000, condition_number(S, Fn, min.val, max.val))
  
  return(mean(sn))
  
}