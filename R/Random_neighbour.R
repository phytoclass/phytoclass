#' Select a random neighbour when the previous random neighbour is beyond 
#' the minimum or maximum value.
#' 
#' @keywords internal
#'
#' @param Fn  xx
#' @param Temp xx
#' @param chlv xx
#' @param s_c xx
#' @param N   xx
#' @param place  xx
#' @param S   xx
#' @param cm   xx
#' @param min.val  xx
#' @param max.val  xx
#'
#' @return
#'
#' @examples
Random_neighbour <- function(Fn, Temp, chlv, s_c, N, place, S, cm, min.val, max.val){
  
  k   <- match(N, place)
  s_c <- s_c[N]
  
  # vectorise function outputs all non-zero elements as a vector
  # (excludes Tchla column)
  wrangled <- Wrangling(Fn, min.val, max.val)
  minF     <- wrangled[[1]][k]
  maxF     <- wrangled[[2]][k]
  SE       <- wrangled[[3]][k]
  
  rand <- round(runif(n = length(s_c), -1, 1), 4)
  SA   <- (SE + Temp * (maxF - minF) * rand)
  
  # F matrix indexes that are outside the range of min and max
  d    <- which(SA < minF | SA > maxF)
  maxd <- maxF
  mind <- minF
  
  loop <- 1
  while (length(d) > 0) {
    loop  <- loop + 1
    nr    <- round(runif(length(d), -1, 1), 4)
    SA2   <- (SE[d] + Temp * maxd[d] - mind[d] * nr)
    SA[d] <- SA2
    d     <- which(SA < minF | SA > maxF)
    
    if (loop > 50) {
      f     <- round(runif(n = length(d), (minF[d] * 1.2), (maxF[d] * 0.80)), 4)
      SA[d] <- f
      d     <- which(SA < minF | SA > maxF)
    }
    
  }
  
  
  Fn           <- Fn[, -ncol(Fn)] #### If error is lower, reassign the values
  Fn[N]        <- SA
  Fn           <- cbind(Fn, chlv)
  colnames(Fn) <- colnames(FALSE)
  
  return(NNLS_MF(Fn, S, cm))
}


