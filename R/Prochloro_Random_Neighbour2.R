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
#' @param chlvp  xx

#'
#' @return
#' 
#' @examples
#' @importFrom stats runif
Prochloro_Random_Neighbour_2 <- function(
  Fn, Temp, chlv, s_c, place, S, cm, min.val, max.val, chlvp
  ) {

  s_c       <- vectorise(s_c[, -ncol(s_c)])
  s_c       <- s_c[-length(s_c)]
  
  wrangle <- Prochloro_Wrangling(Fn, min.val, max.val)
  minF  <- wrangle[[1]]
  maxF  <- wrangle[[2]]
  p_chg <- wrangle[[3]] # pigments to change
  minF  <- minF[-length(minF)]
  maxF  <- maxF[-length(maxF)]
  p_chg <- p_chg[-length(p_chg)]
  
  # randomize pigment ratios
  rand  <- round(runif(n = length(s_c), -1, 1), 4)
  p_new <- p_chg + (Temp) * (maxF - minF) * rand
  oob   <- which(p_new < minF | p_new > maxF)
  
  loop      <- 1
  max_loops <- 100 # Adjust this if needed

  while (length(oob) > 0) {
    loop    <- loop + 1
    nr      <- round(runif(length(oob), -1, 1), 4)
    p_new[oob] <- p_chg[oob] + Temp * (maxF[oob] - minF[oob]) * nr
    oob     <- which(p_new < minF | p_new > maxF)

    # Exit condition after max loops with a safe fallback
    if (loop > max_loops & length(oob) > 0) {
      # sort bound limits for when bounds are small enough to overlap 
      sort_min_max <- cbind(minF[oob] * 1.2, maxF[oob] * 0.80)
      sort_min_max <- t(apply(sort_min_max, 1, sort))
      p_new[oob]   <- round(runif(length(oob), sort_min_max[, 1], sort_min_max[, 2]), 4)
      oob          <- which(p_new < minF | p_new > maxF)
    }
  }

  Fn           <- Fn[, -ncol(Fn)]
  Fn           <- Fn[, -ncol(Fn)]
  Fn[Fn > 0]   <- p_new
  Fn           <- cbind(Fn, chlvp)
  Fn           <- cbind(Fn, chlv)
  colnames(Fn) <- colnames(S)
  return(NNLS_MF(Fn, S, cm))
}
