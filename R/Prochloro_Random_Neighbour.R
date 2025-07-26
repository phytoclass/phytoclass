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
Prochloro_Random_Neighbour <- function(Fn, Temp, chlv, s_c, N, place, S, cm, min.val, max.val)
{
  k <- match(N, place)
  s_c <- s_c[N]
  SE <- Prochloro_Wrangling(Fn, min.val, max.val)[[3]]
  SE <- SE[k]
  minF <- Prochloro_Wrangling(Fn, min.val, max.val)[[1]]
  minF <- minF[k]
  maxF <- Prochloro_Wrangling(Fn, min.val, max.val)[[2]]
  maxF <- maxF[k]
  ki <- maxF - minF
  rand <- round(runif(n = length(s_c), -1, 1), 4)
  SA <- (SE + (Temp) * ki * rand)
  SA <- as.vector(unlist(SA))
  h <- which(maxF == 1)
  if (length(h) > 0) {
    SA[h] <- 1
  }
  d <- which(SA < minF | SA > maxF)
  loop <- 1
  max_loops <- 100
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
    if (loop > max_loops) {
      nn <- (minF[d] + maxF[d])/2
      f <- round(runif(n = length(d), (minF[d] * 1.2),
                       (maxF[d] * 0.8)), 4)
      SA[d] <- f
      d <- which(SA < minF | SA > maxF)
    }
  }
  Fn <- Fn[, 1:ncol(Fn) - 1]
  Fn[N] <- SA
  Fn <- cbind(Fn, chlv)
  colnames(Fn) <- colnames(S)
  F.n <- NNLS_MF(Fn, S, cm)
  return(F.n)
}


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
      # p_new[oob]   <- round(runif(length(oob), minF[oob] * 1.2, maxF[oob] * 0.80), 4)
      oob          <- which(p_new < minF | p_new > maxF)
    }
  }
  Fn           <- Fn[, -ncol(Fn)]
  Fn[N]        <- p_new
  Fn           <- cbind(Fn, chlv)
  colnames(Fn) <- colnames(S)
  return(NNLS_MF(Fn, S, cm))
}

