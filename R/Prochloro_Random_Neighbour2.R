#' Selects a random neighbour for the simulated annealing algorithm.
#' 
#' @keywords internal
#'
#' @param Fn xx
#' @param Temp xx
#' @param chlv xx
#' @param s_c xx
#' @param place  xx
#' @param S  xx
#' @param cm  xx
#' @param min.val  xx
#' @param max.val  xx
#'
#' @return
#' 
#' @examples
#' @importFrom stats runif
Prochloro_Random_Neighbour_2 <- function(Fn, Temp, chlv, s_c, place, S, cm, min.val, max.val) 
  {
  s_c <- phytoclass:::vectorise(s_c[, 1:ncol(s_c) - 1])
  SE <- Prochloro_Wrangling(Fn, min.val, max.val)[[3]]
  minF <- Prochloro_Wrangling(Fn, min.val, max.val)[[1]]
  maxF <- Prochloro_Wrangling(Fn, min.val, max.val)[[2]]
  ki <- maxF - minF
  rand <- round(runif(n = length(s_c), -1, 1), 4)
  SA <- (SE + (Temp) * ki * rand)
  SA <- as.vector(unlist(SA))
  
  d <- which(SA < minF | SA > maxF)
  loop <- 1
  max_loops <- 100 # Adjust this if needed
  
  while (length(d) > 0) {
    loop <- loop + 1
    nr <- round(runif(length(d), -1, 1), 4)
    minr <- as.vector(unlist(minF))
    maxr <- as.vector(unlist(maxF))
    mind <- minr[d]
    maxd <- maxr[d]
    kir <- maxd - mind
    SA2 <- (SE[d] + (Temp) * kir * nr)
    SA[d] <- SA2
    d <- which(SA < minF | SA > maxF)
    
    # Exit condition after max loops with a safe fallback
    if (loop > max_loops) {
      # Set midpoint between minF and maxF as fallback for remaining `SA[d]`
      SA[d] <- (minF[d] + maxF[d]) / 2
      d <- which(SA < minF | SA > maxF) # Recheck bounds after setting fallback
      
      # Break loop if all elements are within bounds
      if (length(d) == 0) { 
        break
      }
    }
  }
  
  Fn <- Fn[, 1:ncol(Fn) - 1]
  Fn[Fn > 0] <- SA
  Fn <- cbind(Fn, chlv)
  colnames(Fn) <- colnames(S)
  F.n <- NNLS_MF(Fn, S, cm)
  return(F.n)
}
