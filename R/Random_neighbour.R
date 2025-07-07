#' Select a random neighbour when the previous random neighbour is beyond 
#' the minimum or maximum value.
#' 
#' @keywords internal
#'
#' @param Fn xx
#' @param Temp xx
#' @param chlv xx
#' @param s_c xx
#' @param N xx
#' @param place xx
#' @param S xx
#' @param S_weights xx
#' @param minF xx
#' @param maxF xx
#'
#' @return
#'
#' @examples
Random_neighbour <- function(Fn, Temp, chlv, s_c, N, place, S, S_weights, minF, maxF) {
  k    <- match(N, place)
  s_c  <- s_c[N]
  SE   <- vectorise(Fn)[k]
  minF <- minF[k]
  maxF <- maxF[k]
  
  rand <- round(runif(n = length(s_c), -1, 1), 4)
  SA   <- SE + (Temp) * (maxF - minF) * rand
  d    <- which(SA < minF | SA > maxF)

  loop <- 1
  while (length(d) > 0) {
    loop  <- loop + 1
    nr    <- round(runif(length(d), -1, 1), 4)
    SA[d] <- SE[d] + Temp * (maxF[d] - minF[d]) * nr
    d     <- which(SA < minF | SA > maxF)
    # print(loop) # if loop is > 3000, just select it from a uniform distribution between the max and min

    if (loop > 50) {
      SA[d] <- round(runif(n = length(d), minF[d] * 1.2, maxF[d] * 0.80), 4)
      d     <- which(SA < minF | SA > maxF)
    }
  }
  
  # If error is lower, reassign the values
  Fn           <- Fn[, -ncol(Fn)] 
  Fn[N]        <- SA
  Fn           <- cbind(Fn, chlv)
  colnames(Fn) <- colnames(FALSE)
  return(NNLS_MF(Fn, S, S_weights))
}
