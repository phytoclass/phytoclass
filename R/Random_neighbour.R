#' Select a random neighbour when the previous random neighbour is beyond 
#' the minimum or maximum value.
#' 
#' @keywords internal
#'
#' @param f_new xx
#' @param Temp xx
#' @param chlv xx
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
Random_neighbour <- function(f_new, Temp, chlv, N, place, S, S_weights, minF, maxF) {
  
  # extract ratios to be changed
  k     <- match(N, place)
  p_chg <- vectorise(f_new)[k] 
  minF  <- minF[k]
  maxF  <- maxF[k]
  
  # randomize pigment ratios
  rand  <- round(runif(n = length(N), -1, 1), 4)
  p_new <- p_chg + (Temp) * (maxF - minF) * rand # new values for ratios
  oob   <- which(p_new < minF | p_new > maxF)    # out of bounds ratios

  loop <- 0
  while (length(oob) > 0) {
    loop  <- loop + 1
    nr    <- round(runif(length(oob), -1, 1), 4)
    p_new[oob] <- p_chg[oob] + Temp * (maxF[oob] - minF[oob]) * nr
    oob   <- which(p_new < minF | p_new > maxF)
    # print(loop) # if loop is > 3000, just select it from a uniform distribution between the max and min

    if (loop > 50 & length(oob) > 0) {
      # sort bound limits for when bounds are small enough to overlap 
      sort_min_max <- cbind(minF[oob] * 1.2, maxF[oob] * 0.80) 
      sort_min_max <- t(apply(sort_min_max, 1, sort))
      p_new[oob]   <- round(runif(length(oob), sort_min_max[, 1], sort_min_max[, 2]), 4)
      oob          <- which(p_new < minF | p_new > maxF)
    }
  }
  
  # If error is lower, reassign the values
  f_new           <- f_new[, -ncol(f_new)] 
  f_new[N]        <- p_new
  f_new           <- cbind(f_new, chlv)
  colnames(f_new) <- colnames(FALSE)
  return(NNLS_MF(f_new, S, S_weights))
}


# ============================================================================ #
# ---- old versions ---- #
# ============================================================================ #

#' #' Selects a random neighbour for the simulated annealing algorithm.
#' #' 
#' #' @keywords internal
#' #'
#' #' @param Fn   xx
#' #' @param Temp xx
#' #' @param chlv xx
#' #' @param s_c xx
#' #' @param place  xx
#' #' @param S      xx
#' #' @param cm   xx
#' #' @param min.val    xx
#' #' @param max.val   xx
#' #'
#' #' @return
#' #'
#' #' @examples
#' #' @importFrom stats runif
#' Random_neighbour2 <- function(Fn, Temp, chlv, s_c, place, S, cm, min.val, max.val){
#'   
#'   s_c  <- vectorise(s_c[, -ncol(s_c)])
#'   
#'   wrangled <- Wrangling(Fn, min.val, max.val)
#'   minF     <- wrangled[[1]]
#'   maxF     <- wrangled[[2]]
#'   SE       <- wrangled[[3]]
#'   
#'   rand <- round(runif(n = length(s_c), -1, 1), 4)
#'   SA   <- SE + Temp * (maxF - minF) * rand
#'   
#'   # F matrix indexes that are outside the range of min and max
#'   d    <- which(SA < minF | SA > maxF) 
#'   
#'   loop <- 1
#'   while (length(d) > 0) {
#'     loop  <- loop + 1
#'     nr    <- round(runif(length(d), -1, 1), 4)
#'     SA[d] <- SE[d] + Temp * (maxF[d] - minF[d]) * nr
#'     d     <- which(SA < minF | SA > maxF)
#'     
#'     if (loop > 50) {
#'       SA[d] <- round(runif(n = length(d), minF[d] * 1.20, maxF[d] * 0.80), 4)
#'       d     <- which(SA < minF | SA > maxF)
#'     }
#'   }
#'   
#'   # If error is lower, reassign the values
#'   Fn           <- Fn[, -ncol(Fn)] 
#'   Fn[Fn > 0]   <- SA
#'   Fn           <- cbind(Fn, chlv)
#'   colnames(Fn) <- colnames(FALSE)
#'   
#'   return(NNLS_MF(Fn, S, cm))
#'   
#' }
