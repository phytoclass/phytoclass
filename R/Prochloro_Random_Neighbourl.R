#' Title
#' 
#' @keywords internal
#'
#' @param Fn 
#' @param Temp 
#' @param chlv 
#' @param s_c 
#' @param place 
#' @param S 
#' @param cm 
#' @param min.val 
#' @param max.val 
#'
#' @return
#' @export
#'
#' @examples
Prochloro_Random_Neighbour_2 <- function (Fn, Temp, chlv, s_c, place, S, cm, min.val, max.val) 
{
  s_c <- phytoclass:::vectorise(s_c[, 1:ncol(s_c) - 1])
  SE <-Prochloro_Wrangling(Fn, min.val, max.val)[[3]]
  minF <- Prochloro_Wrangling(Fn, min.val, max.val)[[1]]
  maxF <- Prochloro_Wrangling(Fn, min.val, max.val)[[2]]
  ki <- maxF - minF
  rand <- round(runif(n = length(s_c), -1, 1), 4)
  SA <- (SE + (Temp) * ki * rand)
  SA <- as.vector(unlist(SA))
  d <- which(SA < minF | SA > maxF)
  length(d)
  loop <- 1
  while (length(d) > 0) {
    loop = loop + 1
    nr <- round(runif(length(d), -1, 1), 4)
    minr <- as.vector(unlist(minF))
    maxr <- as.vector(unlist(maxF))
    mind <- minr[d]
    maxd <- maxr[d]
    kir <- maxd - mind
    SA2 <- (SE[d] + (Temp) * kir * nr)
    SA[d] <- SA2
    d <- which(SA < minF | SA > maxF)
    if (loop > 50) {
      nn <- (minF[d] + maxF[d])/2
      f <- round(runif(n = length(d), (minF[d] * 1.2), 
                       (maxF[d] * 0.8)), 4)
      SA[d] <- f
      d <- which(SA < minF | SA > maxF)
    }
  }
  Fn <- Fn[, 1:ncol(Fn) - 1]
  Fn[Fn > 0] <- SA
  Fn <- cbind(Fn, chlv)
  colnames(Fn) <- colnames(S)
  F.n <- NNLS_MF(Fn, S, cm)
  return(F.n)
}
