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
#' @param cm 
#' @param minF_vec 
#' @param maxF_vec 
#'
#' @return
#'
#' @examples
Prochloro_Random_Neighbour <- function(
    Fn, Temp, chlv, chlvp, s_c, k_idx, S, cm,
    minF_vec, maxF_vec
) {
  core <- Fn[, 1:(ncol(Fn) - 2)]
  core_vec <- as.vector(core)
  
  SE   <- core_vec[k_idx]
  minF <- minF_vec[k_idx]
  maxF <- maxF_vec[k_idx]
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
      SA[d] <- (minF[d] + maxF[d]) / 2
      d <- which(SA < minF | SA > maxF)
    }
  }
  
  v <- core_vec
  v[k_idx] <- SA
  core[] <- v
  
  Fn <- cbind(core, chlvp, chlv)
  colnames(Fn) <- colnames(S)
  
  NNLS_MF(Fn, S, cm)
}
