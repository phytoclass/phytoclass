#' Prochloro random neighbour
#' 
#' @keywords internal
#'
#' @param Fn 
#' @param Temp 
#' @param chlv 
#' @param chlvp
#' @param s_c 
#' @param k_idx 
#' @param S 
#' @param S_weights 
#' @param minF_vec 
#' @param maxF_vec 
#'
#' @return
#'
#' @examples
Prochloro_Random_Neighbour <- function(
    Fn, Temp, chlv, chlvp, s_c, k_idx, S, S_weights,
    minF_vec, maxF_vec
  ) {
  
  core <- Fn[, 1:(ncol(Fn) - 2)]
  core_vec <- as.vector(core)
  
  p_chg <- core_vec[k_idx]
  minF  <- minF_vec[k_idx]
  maxF  <- maxF_vec[k_idx]
  
  # randomize pigment ratios
  rand <- round(runif(n = length(p_chg), -1, 1), 4)
  p_new   <- p_chg + Temp * (maxF - minF) * rand # new values for ratios
  oob <- which(p_new < minF | p_new > maxF)      # out of bounds ratios
  
  loop <- 0
  max_loops <- 100
  while (length(oob) > 0) {
    loop       <- loop + 1
    nr         <- round(runif(length(oob), -1, 1), 4)
    p_new[oob] <- p_chg[oob] + Temp * (maxF[oob] - minF[oob]) * nr
    oob        <- which(p_new < minF | p_new > maxF)
    
    # if (loop > max_loops) {
    #   p_new[oob] <- (minF[oob] + maxF[oob]) / 2
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
  v        <- core_vec
  v[k_idx] <- p_new
  core[]   <- v
  
  f_new           <- cbind(core, chlvp, chlv)
  colnames(f_new) <- colnames(S)
  
  return(NNLS_MF(f_new, S, S_weights))
  
}
