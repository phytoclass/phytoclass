#' Prochloro random neighbour
#' 
#' @keywords internal
#'
#' @param Fn 
#' @param Temp 
#' @param chlv 
#' @param s_c 
#' @param N 
#' @param place 
#' @param S 
#' @param cm 
#' @param min.val 
#' @param max.val 
#'
#' @return
#'
#' @examples
Prochloro_Random_Neighbour <- function(
  Fn, Temp, chlv, s_c, N, place, S, cm, min.val, max.val
  ) {

  k <- match(N, place)
  s_c <- s_c[N]

  wrangle <- Prochloro_Wrangling(Fn, min.val, max.val)
  minF  <- wrangle[[1]][k]
  maxF  <- wrangle[[2]][k]
  p_chg <- wrangle[[3]][k]

  # randomize pigment ratios
  rand <- round(runif(n = length(s_c), -1, 1), 4)
  p_new <- p_chg + Temp * (maxF - minF) * rand

  h <- which(maxF == 1)

  if (length(h) > 0) {
    p_new[h] <- 1
  }

  oob <- which(p_new < minF | p_new > maxF)
  loop <- 1
  max_loops <- 100
  while (length(oob) > 0) {
    loop <- loop + 1
    nr   <- round(runif(length(oob), -1, 1), 4)
    p_new[oob] <- p_chg[oob] + Temp * (maxF[oob] - minF[oob]) * nr
    oob <- which(p_new < minF | p_new > maxF)

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
  Fn[N]        <- p_new
  Fn           <- cbind(Fn, chlv)
  colnames(Fn) <- colnames(S)
  return(NNLS_MF(Fn, S, cm))
}

