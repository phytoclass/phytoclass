#' Calculate the condition number ...
#' 
#' @keywords internal
#'
#' @param S XX     
#' @param Fn  XX   
#' @param min.val XX     
#' @param max.val XX       
#'
#' @return
#'
#' @examples


Condition_test <- function(S, Fn, min.val=NULL, max.val=NULL){
   if (is.null(min.val) & is.null(max.val)) {
    place <- which(Fn > 0)
    K <- Default_min_max(phytoclass::min_max, Fn, place)
    min.val <- K[[1]]
    max.val <- K[[2]]
  }
  condition_number <- function(S, Fd, min.val, max.val) {
    Fz <- vectorise(Fd)
    rand <- vector()
    for (i in 1:length(Fz)) {
      rand[[length(rand) + 1]] <- stats::runif(1, 
                                               min = min.val[i], 
                                               max = max.val[i])
    }
    Fn <- Fd
    Fn[Fn > 0] <- rand
    return(kappa(Fn %*% t(S)))
  }
  sn <- replicate(n = 1000, condition_number(S, Fn, min.val, max.val))
  return(mean(sn))
}
