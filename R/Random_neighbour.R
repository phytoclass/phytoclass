#' Select a random neighbour when the previous random neighbour is beyond 
#' the minimum or maximum value. Part of the simulated annealing algorithm.
#' 
#' @keywords internal
#'
#' @param f_new Current F matrix of pigment ratios
#' @param Temp Current temperature in the annealing process
#' @param chlv Chlorophyll-a column to append to the matrix
#' @param N Indices of pigment ratios to be modified
#' @param place Vector of indices where values are non-zero in the F matrix
#' @param S Matrix of samples (rows) and pigments (columns)
#' @param S_weights Vector of weights for each pigment in NNLS
#' @param minF Minimum bounds for each pigment ratio
#' @param maxF Maximum bounds for each pigment ratio
#'
#' @return A list containing:
#'   \item{F matrix}{The new F matrix with randomly modified values}
#'   \item{RMSE}{Root mean square error of the new solution}
#'   \item{C matrix}{The concentration matrix from NNLS}
#'
#' @examples
#' # Create sample data
#' f_new <- matrix(c(0.5, 0, 0.3, 1,
#'                   0, 0.4, 0.2, 1), 
#'                 nrow=2, byrow=TRUE)
#' Temp <- 0.5
#' chlv <- c(1, 1)
#' N <- c(1, 3)  # indices to modify
#' place <- c(1, 3, 4, 6, 7, 8)  # non-zero indices
#' S <- matrix(runif(12), nrow=4)
#' S_weights <- rep(1, 3)
#' minF <- c(0.3, 0.1)  # for positions in N
#' maxF <- c(0.7, 0.5)  # for positions in N
#' 
#' result <- Random_neighbour(f_new, Temp, chlv, N, place, 
#'                           S, S_weights, minF, maxF)
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

#' Selects a random neighbour for a subset of non-zero pigments that are outside
#' the min and max bounds for the simulated annealing algorithm, specifically 
#' handling Prochlorococcus pigments.
#' 
#' @keywords internal
#'
#' @param f_new F matrix of pigment ratios
#' @param Temp Temperature of the annealing process
#' @param chlv Chlorophyll-a column to append
#' @param chlvp Divinyl chlorophyll-a column for Prochlorococcus
#' @param N Indices of pigment ratios to be modified
#' @param place Vector of indices where values are non-zero in F matrix
#' @param S Matrix of samples (rows) and pigments (columns) 
#' @param S_weights Vector of weights for each pigment in NNLS
#' @param minF Minimum bounds for each pigment ratio
#' @param maxF Maximum bounds for each pigment ratio
#'
#' @return A list containing:
#'   \item{F matrix}{The new F matrix with randomly modified values}
#'   \item{RMSE}{Root mean square error of the new solution}
#'   \item{C matrix}{The concentration matrix from NNLS}
#'
#' @examples
#' # Create sample data with Prochlorococcus
#' f_new <- matrix(c(0.5, 0, 0.3, 0.2, 1,
#'                   0, 0.4, 0.2, 0.3, 1), 
#'                nrow=2, byrow=TRUE)
#' Temp <- 0.5
#' chlv <- c(1, 1)
#' chlvp <- c(0.3, 0.2)  # divinyl chl-a
#' N <- c(1, 3)  # indices to modify
#' place <- c(1, 3, 4, 6, 7, 8, 9, 10)
#' S <- matrix(runif(15), nrow=3)
#' S_weights <- rep(1, 5)
#' minF <- c(0.3, 0.1)  # for positions in N
#' maxF <- c(0.7, 0.5)  # for positions in N
#' 
#' result <- Prochloro_Random_Neighbour(f_new, Temp, chlv, chlvp, N, 
#'                                     place, S, S_weights, minF, maxF)
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
  f_new    <- f_new[, seq(ncol(f_new) - 2)]
  f_new[N] <- p_new
  f_new    <- cbind(f_new, chlvp, chlv)
  
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
#' 
#' #' Selects a random neighbours for all non-zero pigments in the simulated 
#' #' annealing algorithm.
#' #' 
#' #' @keywords internal
#' #'
#' #' @param f_new F matrix
#' #' @param Temp Temperature of the annealing
#' #' @param chlv Chlorophyll-a column
#' #' @param chlvp  Dvchla column
#' #' @param N Indexs of pigment ratios to be changed
#' #' @param place Indexes in F matrix where values are non-zero
#' #' @param S S matrix of samples
#' #' @param S_weights Weights for NNLS algorithm
#' #' @param minF Minimum bounds for each phytoplankton group and pigments
#' #' @param maxF Maximum bounds for each phytoplankton group and pigments
#' #'
#' #' @return
#' #' 
#' #' @examples
#' #' @importFrom stats runif
#' Prochloro_Random_Neighbour_2 <- function(
#'     f_new, Temp, chlv, chlvp, N, place, S, S_weights,
#'     minF, maxF
#' ) {
#'   
#'   # extract ratios to be changed
#'   k     <- match(N, place)
#'   p_chg <- vectorise(f_new)[k] 
#'   minF  <- minF[k]
#'   maxF  <- maxF[k]
#'   
#'   # randomize pigment ratios
#'   rand  <- round(runif(n = length(p_chg), -1, 1), 4)
#'   p_new <- p_chg + Temp * (maxF - minF) * rand # new values for ratios
#'   oob   <- which(p_new < minF | p_new > maxF)  # out of bounds ratios
#'   
#'   loop      <- 0
#'   max_loops <- 100
#'   while (length(oob) > 0) {
#'     loop       <- loop + 1
#'     nr         <- round(runif(length(oob), -1, 1), 4)
#'     p_new[oob] <- p_chg[oob] + Temp * (maxF[oob] - minF[oob]) * nr
#'     oob        <- which(p_new < minF | p_new > maxF)
#'     
#'     # if (loop > max_loops) {
#'     #   p_new[oob] <- (minF[oob] + maxF[oob]) / 2  # midpoint fallback
#'     #   oob <- which(p_new < minF | p_new > maxF)
#'     # }
#'     
#'     if (loop > max_loops & length(oob) > 0) {
#'       # print("on rand neigh loop > 100")
#'       # sort bound limits for when bounds are small enough to overlap 
#'       sort_min_max <- cbind(minF[oob] * 1.2, maxF[oob] * 0.80) 
#'       sort_min_max <- t(apply(sort_min_max, 1, sort))
#'       p_new[oob]   <- round(runif(length(oob), sort_min_max[, 1], sort_min_max[, 2]), 4)
#'       oob          <- which(p_new < minF | p_new > maxF)
#'     }
#'     
#'   }
#'   
#'   # If error is lower, reassign the values
#'   f_new    <- f_new[, seq(ncol(f_new) - 2)]
#'   f_new[N] <- p_new
#'   f_new    <- cbind(f_new, chlvp, chlv)
#'   
#'   return(NNLS_MF(f_new, S, S_weights))
#' }
