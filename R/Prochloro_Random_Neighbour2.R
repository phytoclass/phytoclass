#' Selects a random neighbour for the simulated annealing algorithm.
#' 
#' @keywords internal
#'
#' @param Fn   xx
#' @param Temp xx
#' @param chlv xx
#' @param s_c xx
#' @param place  xx
#' @param S      xx
#' @param cm   xx
#' @param min.val    xx
#' @param max.val   xx
#'
#' @return
#'
#' @examples
#' @importFrom stats runif
Random_neighbour2 <- function(Fn, Temp, chlv, s_c, place, S, cm, min.val, max.val){
  
  s_c  <- vectorise(s_c[, -ncol(s_c)])
  
  wrangled <- Wrangling(Fn, min.val, max.val)
  minF     <- wrangled[[1]]
  maxF     <- wrangled[[2]]
  SE       <- wrangled[[3]]
  
  rand <- round(runif(n = length(s_c), -1, 1), 4)
  SA   <- SE + Temp * (maxF - minF) * rand
  
  # F matrix indexes that are outside the range of min and max
  d    <- which(SA < minF | SA > maxF) 
  
  loop <- 1
  while (length(d) > 0) {
    loop  <- loop + 1
    nr    <- round(runif(length(d), -1, 1), 4)
    SA[d] <- SE[d] + Temp * (maxF[d] - minF[d]) * nr
    d     <- which(SA < minF | SA > maxF)
    
    if (loop > 50) {
      SA[d] <- round(runif(n = length(d), minF[d] * 1.20, maxF[d] * 0.80), 4)
      d     <- which(SA < minF | SA > maxF)
    }
  }
  
  # If error is lower, reassign the values
  Fn           <- Fn[, -ncol(Fn)] 
  Fn[Fn > 0]   <- SA
  Fn           <- cbind(Fn, chlv)
  colnames(Fn) <- colnames(FALSE)
  
  return(NNLS_MF(Fn, S, cm))
  
}
