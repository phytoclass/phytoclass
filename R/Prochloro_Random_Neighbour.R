#' Prochloro random neighbour
#' 
#' @keywords internal
#'
#' @param Fn 
#' @param Temp 
#' @param chlv 
#' @param chlvp
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
    f_new, Temp, chlv, chlvp, N, place, S, S_weights,
    minF, maxF
  ) {
  
  # extract ratios to be changed
  k     <- match(N, place)
  p_chg <- vectorise(f_new)[k] 
  minF  <- minF[k]
  maxF  <- maxF[k]
  
  # randomize pigment ratios
  rand  <- round(runif(n = length(p_chg), -1, 1), 4)
  p_new <- p_chg + Temp * (maxF - minF) * rand # new values for ratios
  oob   <- which(p_new < minF | p_new > maxF)      # out of bounds ratios
  
  loop <- 0
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

  # If error is lower, reassign the values
  f_new           <- f_new[, seq(ncol(f_new) - 2)]
  f_new[N]        <- p_new
  f_new           <- cbind(f_new, chlvp, chlv)

  return(NNLS_MF(f_new, S, S_weights))
  
}
