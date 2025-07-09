#' This is the main phytoclass algorithm. It performs simulated annealing algorithm for S and F matrices. See the examples (Fm, Sm) for how to set up matrices, and the vignette for more detailed instructions. Different pigments and phytoplankton groups may be used.
#' @importFrom progress progress_bar

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
#'
#' @examples
#' # Using the built-in matrices Sm and Fm
#' set.seed(5326)
#' sa.example <- simulated_annealing(Sm, Fm, niter = 5)
#' sa.example$Figure
simulated_annealing <- function(
  S,
  Fmat                 = NULL, 
  user_defined_min_max = NULL,
  do_matrix_checks     = TRUE,
  niter                = 500,
  step                 = 0.009,
  weight.upper.bound   = 30, 
  verbose              = TRUE,
  seed                 = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (is.null(Fmat)) {
    Fmat <- phytoclass::Fm
  }
  
  if (!is.matrix(Fmat)){
    Fmat <- as.matrix(Fmat)
  }
    
  if (is.data.frame(S)) {
    char_cols <- sapply(S, is.character)
    S <- S[, !char_cols]
  }
  
  if (!is.matrix(S)) {
    S <- as.matrix(S)
  }

  if (do_matrix_checks) {
    L    <- Matrix_checks(as.matrix(S), as.matrix(Fmat))
    S    <- as.matrix(L[[1]])
    Fmat <- as.matrix(L[[2]])
  }
  
  S_Chl <- S[, ncol(S)]
  S     <- Normalise_S(S)
  S_weights <- Bounded_weights(S, weight.upper.bound)
  place <- which(Fmat[, -ncol(Fmat)] > 0)
  
  if (is.null(user_defined_min_max)) {
    K <- Default_min_max(phytoclass::min_max, Fmat[, -ncol(Fmat)], place)
  } else {
    K <- Default_min_max(user_defined_min_max,Fmat[, -ncol(Fmat)], place)
  }
    # if (length(min.val) != length(place)) {
    #   message(paste0("\nNumber of rows for user_defined_min_max = ", 
    #                  length(min.val)))
    #   message(paste0("Length of place = ", length(place)))
    #   stop("\nThese do not match.")
    # }
  
  min.val <- K[[1]]
  max.val <- K[[2]]
  
  # ---- start kappa condition check ---- #

  condition.test <- Condition_test(
    S[, -ncol(S)], 
    Fmat[, -ncol(Fmat)], 
    min.val, max.val
    )
  
  if (verbose) {
   message(paste0(
     "\nCondition number = ", round(condition.test),
     "\n\n"
   ))
 }
  
  if (condition.test > 10^5) {
    print("Abort process: condition number of S matrix greater than 100 000\n")
  }
  
  
  # ---- start iteration process ---- #

  Fmat         <- ifelse(Fmat > 0, 1, 0)
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
  
  # extract min and max F matrix thresholds and last column ratio
  wrangled <- Wrangling(f_c, min.val, max.val)
  minF     <- wrangled[[1]]
  maxF     <- wrangled[[2]]
  chlv     <- wrangled[[4]]
  step     <- 1 - step
  
  
  for (k in 1:niter) {
    
    if (!verbose) pb$tick()
    
    Temp <- step^(k)
    
    # needs to be run but not used, due to random number generator
    Random_neighbour(f_c, Temp, chlv, N = place, place, S, S_weights, minF, maxF)
    
    num_loop <- ifelse(k > niter - 20, 300, 120)
    Dn       <- D <- vector("list", num_loop)
    
    for (i in seq(num_loop)) {
      temp_rand <- Random_neighbour(
        f_c, Temp, N = place, chlv, place, S, S_weights, minF, maxF
        )
      D[[i]]  <- temp_rand
      Dn[[i]] <- temp_rand[[2]] # extract RMSE
    }
  
    # select neighbor with lowest RMSE
    nk <- which.min(Dn)
    new_neighbour <- D[[nk]]
    
    num_loop2     <- ifelse(Temp > 0.3, 10, 2)
    new_neighbour <- Steepest_Descent(new_neighbour[[1]], place, S, S_weights, 
                                      num.loops = num_loop2)
    
    f_n      <- new_neighbour[[1]]

    # check if ratios are out of bounds (min\max)
    oob <- which(vectorise(f_n[, -ncol(f_n)]) < minF | 
                 vectorise(f_n[, -ncol(f_n)]) > maxF)
    
    # if new lowest neighbor has a ratio outside of bounds, will change only 
    # those ones
    while (length(oob) > 0) {
      
      N <- place[oob] # where in F matrix is the ratio out of bounds
      Dn2 <- D2 <- vector("list", num_loop)

      for (i in seq(num_loop)) {
        temp_rand <- Random_neighbour(
          f_n, Temp, chlv, N, place, S, S_weights,minF, maxF
          )
        D2[[i]]  <- temp_rand
        Dn2[[i]] <- temp_rand[[2]] # extract RMSE
      }
      
      # select new neighbor with lowest RMSE 
      nk <- which.min(c(Dn, Dn2))
      
      new_neighbour <- c(D, D2)[[nk]]
      f_n           <- new_neighbour[[1]]
      
      # check if ratios are out of bounds (min\max)
      oob <- which(vectorise(f_n[, -ncol(f_n)]) < minF |
                   vectorise(f_n[, -ncol(f_n)]) > maxF) 
      
    } 
    
    # check RMSE of neighbor is better than current 
    f_n_err  <- new_neighbour[[2]] # RMSE
    if (f_n_err < f_c_err || 
        exp(-(f_n_err - f_c_err)) < stats::runif(1, 0, 1)
        ) {
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
          "Current error: ", round(f_c_err, 4),
          "\nNeighbour's error: ", round(f_n_err, 4),
          "\nTemperature (%): ", round(Temp * 100, 2), "\n"
          )
        )
    }
  }

  final_results <- NNLS_MF_Final(f_b, S, S_Chl, S_weights)
  return(final_results)
}
