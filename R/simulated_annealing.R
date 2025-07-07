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
  cm    <- Bounded_weights(S, weight.upper.bound)
  place <- which(Fmat[,1:ncol(Fmat) - 1] > 0)
  
  if (is.null(user_defined_min_max)) {
    K <- Default_min_max(phytoclass::min_max, Fmat[,1:ncol(Fmat) - 1], place)
    min.val <- K[[1]]
    max.val <- K[[2]]
  } else {
    K <- Default_min_max(user_defined_min_max,Fmat[, 1:ncol(Fmat) - 1], place)
    min.val <- K[[1]]
    max.val <- K[[2]]
    # if (length(min.val) != length(place)) {
    #   message(paste0("\nNumber of rows for user_defined_min_max = ", 
    #                  length(min.val)))
    #   message(paste0("Length of place = ", length(place)))
    #   stop("\nThese do not match.")
    # }
  }
  
  # start kappa condition check
  condition.test <- Condition_test(
    S[,1:ncol(S) - 1], 
    Fmat[,1:ncol(Fmat) - 1], 
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

  Fmat <- ifelse(Fmat > 0, 1, 0)
  # SE   <- vectorise(Fmat)
  nnls_initial <- NNLS_MF(Fmat, S, cm)
  
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
  
  chlv    <- Wrangling(f_c, min.val, max.val)[[4]]
  
  for (k in 1:niter) {
    
    if (!verbose) pb$tick()
    
    Temp <- (1 - step)^(k)
    # chlv <- Wrangling(f_c, min.val, max.val)[[4]]
    new_neighbour <- Random_neighbour(f_c, Temp, chlv, f_c, N = place,
                                       place, S, cm, min.val, max.val)
    
    num_loop <- ifelse(k > niter - 20, 300, 120)
    Dn       <- D <- vector("list", num_loop)
    
    for (i in seq(num_loop)) {
      # chlv    <- Wrangling(f_c, min.val, max.val)[[4]]
      temp    <- Random_neighbour(f_c, Temp, N = place, chlv, f_c, place, S, cm, min.val, max.val)
      D[[i]]  <- temp
      Dn[[i]] <- temp[[2]] # extract RMSE
    }

    nk <- which.min(Dn)
    new_neighbour <- D[[nk]]
    
    num_loop2     <- ifelse(Temp > 0.3, 10, 2)
    new_neighbour <- SAALS(new_neighbour[[1]], min.val, 
                           max.val, place, S, cm, num.loops = num_loop2)
    
    wrangled <- Wrangling(new_neighbour[[1]], min.val, max.val)
    minF     <- wrangled[[1]]
    maxF     <- wrangled[[2]]
    
    f_n      <- new_neighbour[[1]]
    f_n_err  <- new_neighbour[[2]]
    loop     <- 1
    d        <- which(vectorise(f_n[, 1:(ncol(f_n) - 1)]) < minF | 
                      vectorise(f_n[, 1:(ncol(f_n) - 1)]) > maxF)
    
    while (length(d) > 0) {
      
      num_loop3 <- ifelse(k > niter - 20, 300, 120)
      N <- place[d]
      
      Dn2 <- D2 <- vector("list", num_loop3)
      for (i in seq(num_loop3)) {
        # chlv     <- Wrangling(f_n, min.val, max.val)[[4]]
        temp     <- Random_neighbour(f_n, Temp, chlv, f_n, N, place, S, cm, min.val, max.val)
        D2[[i]]  <- temp
        Dn2[[i]] <- temp[[2]] # extract RMSE
      }
      
      nk <- which.min(c(Dn, Dn2))
      
      new_neighbour <- c(D, D2)[[nk]]
      f_n           <- new_neighbour[[1]]
      f_n_err       <- new_neighbour[[2]]
      
      wrangled <- Wrangling(new_neighbour[[1]], min.val, max.val)
      minF     <- wrangled[[1]]
      maxF     <- wrangled[[2]]
      
      d <- which(vectorise(f_n[,1:(ncol(f_n) - 1)]) < minF |
                 vectorise(f_n[,1:(ncol(f_n) - 1)]) > maxF) 
      
    } 
    
    # check RMSE
    # A <-  target(f_n_err)/target(f_c_err)
    diff <- f_n_err - f_c_err
    if (f_n_err < f_c_err || 
        exp(-(f_n_err - f_c_err)) < stats::runif(1, 0, 1)
        ) {
      f_c     <- f_n
      f_c_err <- f_n_err
    }
    
    if (verbose) {
      message(paste("Current error: ", round(f_c_err, 4)))
      message(paste("Neighbour's error: ", round(f_n_err, 4)))
      message(paste("Temperature (%): ", round(Temp * 100, 2)))
      message(" ")
    }
    
    if (f_n_err < f_b_err) {
      f_b     <- f_n
      f_b_err <- f_n_err
    }

  }
  
  res <- list(f_b, f_b_err)
  A   <- res[[1]]
  
  final.results <- NNLS_MF_Final(A, S, S_Chl, cm)
  return(final.results)
}
