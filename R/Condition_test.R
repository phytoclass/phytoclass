#' Calculate the condition number ...
#'
#' @param Fn
#' @param min.val
#' @param max.val
#'
#' @return
#' @export
#'
#' @examples


Condition_test <- function(Fn, min.val, max.val){
  
  condition_number <- function(Fd, min.val, max.val){
    Fz <- vectorise(Fd)
    rand <- vector()
    for (i in 1:length(Fz)){
      rand[[length(rand)+1]] <- stats::runif(1, min=min.val[i], max=max.val[i])
    }
    Fn <- Fd
    Fn[Fn >0] <- rand
    return(kappa(Fn %*% t(S)))
  }
  
  sn <- replicate(n = 10000, condition_number(Fn, min.val, max.val))
  
  return(mean(sn))
  
}