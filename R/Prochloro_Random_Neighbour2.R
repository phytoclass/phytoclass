#' Selects a random neighbour for the simulated annealing algorithm.
#' 
#' @keywords internal
#'
#' @param Fn xx
#' @param Temp xx
#' @param chlv xx
#' @param S  xx
#' @param S_weights  xx
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
    Fn, Temp, chlv, S, S_weights,
    minF_vec, maxF_vec, chlvp, Fi_mask
  ) {
  core <- Fn[, seq(ncol(Fn) - 2)]
  mask_core <- Fi_mask[, seq(ncol(Fi_mask) - 2)]
  
  core_vec <- as.vector(core)
  mask_vec <- as.vector(mask_core)
  idx      <- which(mask_vec > 0)
  
  p_chg <- core_vec[idx]
  minF  <- minF_vec[idx]
  maxF  <- maxF_vec[idx]
  
  # randomize pigment ratios
  rand  <- round(runif(n = length(p_chg), -1, 1), 4)
  p_new <- p_chg + Temp * (maxF - minF) * rand # new values for ratios
  oob   <- which(p_new < minF | p_new > maxF)  # out of bounds ratios
  
  loop      <- 0
  max_loops <- 100
  while (length(oob) > 0) {
    loop       <- loop + 1
    nr         <- round(runif(length(oob), -1, 1), 4)
    p_new[oob] <- p_chg[oob] + Temp * (maxF[oob] - minF[oob]) * nr
    oob        <- which(p_new < minF | p_new > maxF)
    
    # if (loop > max_loops) {
    #   p_new[oob] <- (minF[oob] + maxF[oob]) / 2  # midpoint fallback
    #   oob <- which(p_new < minF | p_new > maxF)
    # }
    
    if (loop > max_loops & length(oob) > 0) {
      # sort bound limits for when bounds are small enough to overlap 
      sort_min_max <- cbind(minF[oob] * 1.2, maxF[oob] * 0.80) 
      sort_min_max <- t(apply(sort_min_max, 1, sort))
      p_new[oob]   <- round(runif(length(oob), sort_min_max[, 1], sort_min_max[, 2]), 4)
      oob          <- which(p_new < minF | p_new > maxF)
    }
    
  }
  
  # Write back ONLY to allowed cells
  v      <- core_vec
  v[idx] <- p_new
  core[] <- v
  
  F_new <- cbind(core, chlvp, chlv)
  colnames(F_new) <- colnames(S)
  
  return(NNLS_MF(F_new, S, S_weights))
}
