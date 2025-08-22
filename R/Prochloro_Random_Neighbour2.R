#' Selects a random neighbour for the simulated annealing algorithm.
#' 
#' @keywords internal
#'
#' @param Fn xx
#' @param Temp xx
#' @param chlv xx
#' @param s_c xx
#' @param place_full  xx
#' @param S  xx
#' @param cm  xx
#' @param minF_vec  xx
#' @param maxF_vec  xx
#' @param chlvp  xx
#' @param Fi_mask  xx

#'
#' @return
#' 
#' @examples
#' @importFrom stats runif
Prochloro_Random_Neighbour_2 <- function(
    Fn, Temp, chlv, s_c, place_full, S, cm,
    minF_vec, maxF_vec, chlvp, Fi_mask
) {
  core <- Fn[, 1:(ncol(Fn) - 2)]
  mask_core <- Fi_mask[, 1:(ncol(Fi_mask) - 2)]
  
  core_vec <- as.vector(core)
  mask_vec <- as.vector(mask_core)
  idx <- which(mask_vec > 0)
  
  SE   <- core_vec[idx]
  minF <- minF_vec[idx]
  maxF <- maxF_vec[idx]
  ki   <- maxF - minF
  
  rand <- round(runif(n = length(SE), -1, 1), 4)
  SA   <- SE + Temp * ki * rand
  
  d <- which(SA < minF | SA > maxF)
  loop <- 1; max_loops <- 100
  while (length(d) > 0) {
    loop <- loop + 1
    nr <- round(runif(length(d), -1, 1), 4)
    SA[d] <- SE[d] + Temp * (maxF[d] - minF[d]) * nr
    d <- which(SA < minF | SA > maxF)
    if (loop > max_loops) {
      SA[d] <- (minF[d] + maxF[d]) / 2  # midpoint fallback
      d <- which(SA < minF | SA > maxF)
    }
  }
  
  # Write back ONLY to allowed cells
  v <- core_vec
  v[idx] <- SA
  core[] <- v
  
  Fn <- cbind(core, chlvp, chlv)
  colnames(Fn) <- colnames(S)
  
  NNLS_MF(Fn, S, cm)
}
