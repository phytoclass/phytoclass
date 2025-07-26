#' Perform simulated annealing algorithm for samples with divinyl chlorophyll and prochlorococcus. Divinyl chlorophyll must be the final column of both S and F matrices, with chlorophyll a the 2nd to last column. See how the example Sp and Fp matrices are organised.
#'
#' @param S   Sample data matrix – a matrix of pigment samples
#' @param Fmat   Pigment to Chl a matrix
#' @param user_defined_min_max data frame with some format as min_max built-in data
#' @param do_matrix_checks     This should only be set to TRUE when using the default values. This will remove pigment columns that have column sums of 0. Set to FALSE if using customised names for pigments and phytoplankton groups
#' @param niter Number of iterations (default is 500)
#' @param step  Step ratio used (default is 0.009)
#' @param weight.upper.bound Upper limit of the weights applied (default value is 30). 
#' @param verbose Logical value. Output error and temperature at each iteration. Default value of TRUE
#'
#' @return A list containing 
#' \enumerate{
#'  \item Fmat matrix
#'  \item RMSE (Root Mean Square Error)
#'  \item condition number
#'  \item Class abundances
#'  \item Figure (plot of results)
#'  \item MAE (Mean Absolute Error)
#'  \item Error
#'  }
#' @export
#'
#' @examples
#' # Using the built-in matrices Sp and Fp.
#' set.seed(5326)
#' sa.example <- simulated_annealing_Prochloro(Sp, Fp, niter = 1)
#' sa.example$Figure
# simulated_annealing_Prochloro <- function (S, 
#                                            Fmat = NULL,
#                                            user_defined_min_max = NULL,
#                                            do_matrix_checks = TRUE, 
#                                            niter = 500, 
#                                            step = 0.009, 
#                                            weight.upper.bound = 30, 
#                                            verbose = TRUE) 
# {
#   if (is.null(Fmat)) {
#     Fmat <- phytoclass::Fp
#   }
#   if (is.data.frame(S)) {
#     char_cols <- sapply(S, is.character)
#     S <- S[, !char_cols]
#   }
#   if (do_matrix_checks) {
#     L <- Matrix_checks(as.matrix(S), as.matrix(Fmat))
#     S <- as.matrix(L[[1]])
#     Fmat <- as.matrix(L[[2]])
#   }
#   S_Chl <- S[, ncol(S)]
#   S_dvChl <- S[, ncol(S) - 1]
#   S <- Normalise_S(S)
#   cm <- Bounded_weights(S, weight.upper.bound)
#   place <- which(Fmat[, 1:ncol(Fmat) - 1] > 0)
#   if (is.null(user_defined_min_max)) {
#     K <- Default_min_max(phytoclass::min_max, 
#                                       Fmat[, 1:ncol(Fmat) - 1], place)
#     min.val <- K[[1]]
#     max.val <- K[[2]]
#   }
#   else {
#     K <- Default_min_max(user_defined_min_max, 
#                                       Fmat[, 1:ncol(Fmat) - 1], place)
#     min.val <- K[[1]]
#     max.val <- K[[2]]
#   }
#   condition.test <- Condition_test(S[, 1:ncol(S) - 
#                                                     1], Fmat[, 1:ncol(Fmat) - 1], min.val, max.val)
#   if (verbose) {
#     message(paste0("\nCondition number = ", round(condition.test), 
#                    "\n\n"))
#   }
#   if (condition.test > 10^5) {
#     print("Abort process: condition number of S matrix greater than 100 000\n")
#   }
#   Fi <- ifelse(Fmat > 0, 1, 0)
#   SE <- vectorise(Fi)
#   nc <- NNLS_MF(Fi, S, cm)
#   s_b <- s_c <- s_n <- nc[[1]]
#   f_b <- f_c <- f_n <- nc[[2]]
#   for (k in 1:niter) {
#     Temp <- (1 - step)^(k)
#     chlv <- Prochloro_Wrangling(s_c, min.val, max.val)[[4]]
#     chlvp <- Prochloro_Wrangling(s_c, min.val, max.val)[[5]]
#     new_neighbour <- Prochloro_Random_Neighbour_2(s_c, Temp,
#                                                   chlv, s_c, place, S, cm, min.val, max.val, chlvp)
#     # new_neighbour <- Prochloro_Random_Neighbour(s_c, Temp, 
#     #                                             chlv, s_c, N = place, place, 
#     #                                             S, cm, min.val, max.val, chlvp)
#     D <- list()
#     if (k > niter - 20) {
#       for (i in 1:300) {
#         chlv <- Prochloro_Wrangling(s_c, min.val, max.val)[[4]]
#         chlvp <- Prochloro_Wrangling(s_c, min.val, max.val)[[5]]
#         D[[length(D) + 1]] <- Prochloro_Random_Neighbour_2(s_c,
#                                                            Temp, chlv, s_c, place, S, cm, min.val, max.val,
#                                                            chlvp)
#         # D[[length(D) + 1]] <- Prochloro_Random_Neighbour(s_c, Temp, chlv, s_c, 
#         #                                                    N = place, place, S, 
#         #                                                    cm, min.val, max.val, 
#         #                                                    chlvp)
#       }
#     }
#     else {
#       for (i in 1:120) {
#         chlv <- Prochloro_Wrangling(s_c, min.val, max.val)[[4]]
#         chlvp <- Prochloro_Wrangling(s_c, min.val, max.val)[[5]]
#         D[[length(D) + 1]] <- Prochloro_Random_Neighbour_2(s_c,
#                                                            Temp, chlv, s_c, place, S, cm, min.val, max.val,
#                                                            chlvp)
#         # D[[length(D) + 1]] <- Prochloro_Random_Neighbour(s_c, Temp, chlv, s_c, 
#         #                                                  N = place, place, S, cm, 
#         #                                                  min.val, max.val, chlvp)
#       }
#     }
#     Dn <- list()
#     for (i in D) {
#       Dn[[length(Dn) + 1]] <- i[[2]]
#     }
#     nk <- which.min(Dn)
#     new_neighbour <- D[[nk]]
#     if (Temp > 0.3) {
#       new_neighbour <- SAALS(new_neighbour[[1]], 
#                                           min.val, max.val, place, S, cm, num.loops = 10)
#     }
#     else {
#       new_neighbour <- SAALS(new_neighbour[[1]], 
#                                           min.val, max.val, place, S, cm, num.loops = 2)
#     }
#     minF <- Prochloro_Wrangling(new_neighbour[[1]], min.val, 
#                                 max.val)[[1]]
#     maxF <- Prochloro_Wrangling(new_neighbour[[1]], min.val, 
#                                 max.val)[[2]]
#     s_n <- new_neighbour[[1]]
#     f_n <- new_neighbour[[2]]
#     loop <- 1
#     d <- which(vectorise(s_n[, 1:(ncol(s_n) - 
#                                                  2)]) < minF[1:length(minF) - 1] | vectorise(s_n[, 
#                                                                                                               1:(ncol(s_n) - 2)]) > maxF[1:length(maxF) - 1])
#     while (length(d) > 0) {
#       if (k > niter - 20) {
#         N <- place[d]
#         for (i in 1:300) {
#           chlv <- Prochloro_Wrangling(s_n, min.val, max.val)[[4]]
#           chlvp <- Prochloro_Wrangling(s_c, min.val, 
#                                        max.val)[[5]]
#           D[[length(D) + 1]] <- Prochloro_Random_Neighbour(s_n, 
#                                                            Temp, chlv, s_n, N, place, S, cm, min.val, 
#                                                            max.val)
#         }
#         Dn <- list()
#         for (i in D) {
#           Dn[[length(Dn) + 1]] <- i[[2]]
#         }
#         nk <- which.min(Dn)
#         new_neighbour <- D[[nk]]
#         s_n <- new_neighbour[[1]]
#         f_n <- new_neighbour[[2]]
#         minF <- Prochloro_Wrangling(new_neighbour[[1]], 
#                                     min.val, max.val)[[1]]
#         maxF <- Prochloro_Wrangling(new_neighbour[[1]], 
#                                     min.val, max.val)[[2]]
#         d <- which(vectorise(s_n[, 1:(ncol(s_n) - 
#                                                      2)]) < minF[1:length(minF) - 1] | vectorise(s_n[, 
#                                                                                                                   1:(ncol(s_n) - 2)]) > maxF[1:length(maxF) - 
#                                                                                                                                                1])
#       }
#       else {
#         N <- place[d]
#         for (i in 1:120) {
#           chlv <- Prochloro_Wrangling(s_n, min.val, max.val)[[4]]
#           chlvp <- Prochloro_Wrangling(s_c, min.val, 
#                                        max.val)[[5]]
#           D[[length(D) + 1]] <- Prochloro_Random_Neighbour(s_n, 
#                                                            Temp, chlv, s_n, N, place, S, cm, min.val, 
#                                                            max.val)
#         }
#         Dn <- list()
#         for (i in D) {
#           Dn[[length(Dn) + 1]] <- i[[2]]
#         }
#         nk <- which.min(Dn)
#         new_neighbour <- D[[nk]]
#         s_n <- new_neighbour[[1]]
#         f_n <- new_neighbour[[2]]
#         minF <- Prochloro_Wrangling(new_neighbour[[1]], 
#                                     min.val, max.val)[[1]]
#         maxF <- Prochloro_Wrangling(new_neighbour[[1]], 
#                                     min.val, max.val)[[2]]
#         d <- which(vectorise(s_n[, 1:(ncol(s_n) - 
#                                                      2)]) < minF[1:length(minF) - 1] | vectorise(s_n[, 
#                                                                                                                   1:(ncol(s_n) - 2)]) > maxF[1:length(maxF) - 
#                                                                                                                                                1])
#       }
#     }
#     diff <- f_n - f_c
#     if (f_n < f_c || exp(-(f_n - f_c)) < stats::runif(1, 
#                                                       0, 1)) {
#       s_c <- s_n
#       f_c <- f_n
#     }
#     if (verbose) {
#       message(paste("Current error: ", round(f_c, 4)))
#       message(paste("Neighbour's error: ", round(f_n, 4)))
#       message(paste("Temperature (%): ", round(Temp * 100, 
#                                                2)))
#       message(" ")
#     }
#     if (f_n < f_b) {
#       s_b <- s_n
#       f_b <- f_n
#     }
#   }
#   res <- list(s_b, f_b)
#   A <- res[[1]]
#   final.results <- Prochloro_NNLS_MF_Final(A, S, S_Chl, cm, 
#                                            S_dvChl)
#   return(final.results)
# }


#-------------------------------------------------------------------------- 
# clean up

simulated_annealing_Prochloro <- function(
    S, 
    Fmat                 = NULL,
    user_defined_min_max = NULL,
    do_matrix_checks     = TRUE, 
    niter                = 500, 
    step                 = 0.009, 
    weight.upper.bound   = 30, 
    verbose              = TRUE,
    seed                 = NULL) {
  
  # set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (is.null(Fmat)) {
    Fmat <- phytoclass::Fp
  }
  
  if (is.data.frame(S)) {
    char_cols <- sapply(S, is.character)
    S         <- S[, !char_cols]
  }
  
  if (!is.matrix(S)) {
    S <- as.matrix(S)
  }
  
  if (do_matrix_checks) {
    mat_n <- Matrix_checks(as.matrix(S), as.matrix(Fmat))
    S     <- as.matrix(mat_n[[1]])
    Fmat  <- as.matrix(mat_n[[2]])
  }
  
  S_Chl     <- S[, ncol(S)] # chla 
  S_dvChl   <- S[, ncol(S) - 1] # dvchla
  S         <- Normalise_S(S)
  S_weights <- Bounded_weights(S, weight.upper.bound)
  place     <- which(Fmat[, -ncol(Fmat)] > 0) # non-zero, non-chla pigments
  
  if (is.null(user_defined_min_max)) {
    min_max_mat <- Default_min_max(phytoclass::min_max, Fmat[, -ncol(Fmat)], place)
  } else {
    min_max_mat <- Default_min_max(user_defined_min_max, Fmat[, -ncol(Fmat)], place)
    }

  # ---- start kappa condition check ---- #
  condition.test <- Condition_test(
    S[, -ncol(S)], 
    Fmat[, -ncol(Fmat)], 
    min_max_mat[[1]], min_max_mat[[2]]
  )
  
  if (verbose) {
    message(paste0("\nCondition number = ", round(condition.test), "\n\n"))
  }
  
  if (condition.test > 10^5) {
    print("Abort process: condition number of S matrix greater than 100 000\n")
  }
  
  
  # ---- start iteration process ---- #
  Fmat <- ifelse(Fmat > 0, 1, 0)
  nnls_initial <- NNLS_MF(Fmat, S, S_weights)
  
  # initialize F matrix and RMSE
  # n = neighbor
  # c = current
  # b = best
  f_b     <- f_c     <- f_n     <- nnls_initial[[1]] # F matrix
  f_b_err <- f_c_err <- f_n_err <- nnls_initial[[2]] # RMSE
  
  
  # Initialize progress bar if verbose is FALSE
  if (!verbose) {
    pb <- progress::progress_bar$new(
      format = "  Simulated Annealing [:bar] :percent ETA: :eta",
      total  = niter, clear = FALSE, width = 60
    )
  }
  
  step <- 1 - step
  
  for (k in 1:niter) {
    
    if (!verbose) pb$tick()
    
    Temp    <- step^k
    wrangle <- Prochloro_Wrangling(f_c, min_max_mat[[1]], min_max_mat[[2]])
    chlv    <- wrangle[[4]]
    chlvp   <- wrangle[[5]]
    
    # needs to be run but not used, due to random number generator
    Prochloro_Random_Neighbour_2(
      f_c, Temp, chlv, f_c, place, S, S_weights, min_max_mat[[1]], min_max_mat[[2]], chlvp
      )
    
    num_loop <- ifelse(k > niter - 20, 300, 120)
    
    # set list to store random neighbors
    rand_itr_err <- rand_itr <- vector("list", num_loop) 

    for (i in seq(num_loop)) {
      wrangle   <- Prochloro_Wrangling(f_c, min_max_mat[[1]], min_max_mat[[2]])
      chlv      <- wrangle[[4]]
      chlvp     <- wrangle[[5]]
      temp_rand <- Prochloro_Random_Neighbour_2(
        f_c, Temp, chlv, f_c, place, S, S_weights, 
        min_max_mat[[1]], min_max_mat[[2]], chlvp
        )
      
      rand_itr[[i]]     <- temp_rand
      rand_itr_err[[i]] <- temp_rand[[2]] # extract RMSE
      
    }

    # select neighbor with lowest RMSE
    low_indx      <- which.min(rand_itr_err)
    new_neighbour <- rand_itr[[low_indx]]
    
    # steepest descent
    num_loop2     <- ifelse(Temp > 0.3, 10, 2)
    new_neighbour <- Steepest_Descent(new_neighbour[[1]], place, S, S_weights, 
                                      num_loop2)
    
    f_n     <- new_neighbour[[1]]
    
    # check if ratios are out of bounds (min\max)
    vect    <- vectorise(f_n[, seq(ncol(f_n) - 2)])
    wrangle <- Prochloro_Wrangling(new_neighbour[[1]], min_max_mat[[1]], min_max_mat[[2]])
    minF    <- wrangle[[1]]
    maxF    <- wrangle[[2]]
    oob     <- which(vect < minF[-length(minF)] | vect > maxF[-length(maxF)])
    
    while (length(oob) > 0) {
      
      oob_indx  <- place[oob] # where in F matrix is the ratio out of bounds
      num_loop  <- ifelse(k > niter - 20, 300, 120)
      
      # set list to store random neighbors
      rand_itr2_err <- rand_itr2 <- vector("list", num_loop)
      
      for (i in seq(num_loop)) {
        chlv      <- Prochloro_Wrangling(f_n, min_max_mat[[1]], min_max_mat[[2]])[[4]]
        chlvp     <- Prochloro_Wrangling(f_c, min_max_mat[[1]], min_max_mat[[2]])[[5]]
        temp_rand <- Prochloro_Random_Neighbour(
          f_n, Temp, chlv, f_n, oob_indx, place, S, S_weights, 
          min_max_mat[[1]], min_max_mat[[2]]
          )
        rand_itr2[[i]]     <- temp_rand
        rand_itr2_err[[i]] <- temp_rand[[2]] # extract RMSE
        
      }
      
      # select new neighbor with lowest RMSE 
      low_indx      <- which.min(c(rand_itr_err, rand_itr2_err))
      new_neighbour <- c(rand_itr, rand_itr2)[[low_indx]]
      
      f_n     <- new_neighbour[[1]]
      
      # check if ratios are out of bounds (min\max)
      vect    <- vectorise(f_n[, seq(ncol(f_n) - 2)])
      wrangle <- Prochloro_Wrangling(new_neighbour[[1]], min_max_mat[[1]], min_max_mat[[2]])
      minF    <- wrangle[[1]]
      maxF    <- wrangle[[2]]
      oob     <- which(vect < minF[-length(minF)] | vect > maxF[-length(maxF)])
      
    }
    
    # check RMSE of neighbor is better than current 
    f_n_err  <- new_neighbour[[2]] # RMSE
    if (f_n_err < f_c_err || exp(-(f_n_err - f_c_err)) < stats::runif(1, 0, 1)) {
      f_c     <- f_n
      f_c_err <- f_n_err
    }
    
    # update F matrix and RMSE when neighbor F matrix has lower RMSE
    if (f_n_err < f_b_err) {
      f_b     <- f_n
      f_b_err <- f_n_err
    }
    
    if (verbose) {
      message(
        paste(
          "Iterations:        ", sprintf("%03d", k), "of", niter,
          "\nCurrent error:     ", round(f_c_err, 4),
          "\nNeighbour's error: ", round(f_n_err, 4),
          "\nTemperature (%):   ", round(Temp * 100, 2), "\n"
        )
      )
    }
  }

  final_results <- Prochloro_NNLS_MF_Final(f_b, S, S_Chl, S_weights, S_dvChl)
  
  return(final_results)
}

