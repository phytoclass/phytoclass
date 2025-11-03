#' This is the main phytoclass algorithm. It performs simulated annealing
#' algorithm for S and F matrices. See the examples (Fm, Sm) for how to set
#' up matrices, and the vignette for more detailed instructions. Different
#' pigments and phytoplankton groups may be used.
#' @importFrom progress progress_bar

#' @param S   Sample data matrix – a matrix of pigment samples
#' @param Fmat   Pigment to Chl a matrix
#' @param user_defined_min_max data frame with some format as min_max built-in data
#' @param do_matrix_checks This should only be set to TRUE when using the default
#'                         values. This will remove pigment columns that have 
#'                         column sums of 0. Set to FALSE if using customised 
#'                         names for pigments and phytoplankton groups
#' @param niter Number of iterations (default is 500)
#' @param step  Step ratio used (default is 0.009)
#' @param weight.upper.bound Upper limit of the weights applied (default value is 30). 
#' @param verbose Logical value. Output error and temperature at each
#'   iteration. Default value of TRUE
#' @param seed Set number to reproduce the same results
#' @param check_converge TRUE/FALSE/integer; set the number of F matrices to 
#'                        for convergence checking
#' @param alt_pro_name Optional: additional alternate versions of 
#'                     divinyl-chlorophyll-a spellings used to detect 
#'                     prochlorococcus (Default: "dvchl", "dvchla", "dv_chla")
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
#'  \item F_mat_iter
#'  \item converge_plot 
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
  seed                 = NULL,
  check_converge       = 100,
  alt_pro_name         = NULL
  ) {

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
    mat_check <- Matrix_checks(as.matrix(S), as.matrix(Fmat))
    S    <- as.matrix(mat_check[[1]])
    Fmat <- as.matrix(mat_check[[2]])
  }
  
  # Check for dvchl/dvchla
  col_names       <- tolower(colnames(S))
  pro_name        <- tolower(c("dvchl", "dvchla", "dv_chla", alt_pro_name))
  penultimate_col <- col_names %in% pro_name
  
  if (any(penultimate_col)) {
    message("Detected dvchl column. Using simulated_annealing_Prochloro().")
    pro_row <- which(Fmat[,colnames(Fmat)[penultimate_col]] > 0)

    if (pro_row != nrow(Fmat)) {
      message("Moving phytoplankton group with Dvchla to last row.")
      Fmat <- Fmat[c(setdiff(seq_len(nrow(Fmat)), pro_row), pro_row), ]
    }
    return(
      simulated_annealing_Prochloro(
        S    = S,
        Fmat = Fmat,
        user_defined_min_max = user_defined_min_max,
        do_matrix_checks     = do_matrix_checks,
        niter = niter,
        step  = step,
        weight.upper.bound = weight.upper.bound,
        verbose = verbose
        )
      )
  }
  
  S_Chl <- S[, ncol(S)]
  S     <- Normalise_S(S)
  S_weights <- Bounded_weights(S, weight.upper.bound)
  place <- which(Fmat[, -ncol(Fmat)] > 0)
  
  if (is.null(user_defined_min_max)) {
    min_max_mat <- Default_min_max(phytoclass::min_max, Fmat[, -ncol(Fmat)])
  } else {
    min_max_mat <- Default_min_max(user_defined_min_max,Fmat[, -ncol(Fmat)])
  }
    # if (length(min_max_mat[[1]]) != length(place)) {
    #   message(paste0("\nNumber of rows for user_defined_min_max = ", 
    #                  length(min_max_mat[[1]])))
    #   message(paste0("Length of place = ", length(place)))
    #   stop("\nThese do not match.")
    # }
  
  # ---- start kappa condition check ---- #

  condition.test <- Condition_test(
    S[, -ncol(S)], 
    Fmat[, -ncol(Fmat)], 
    min_max_mat[[1]], min_max_mat[[2]]
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
  
  # initialize progress bar if verbose is FALSE
  if (!verbose) {
    pb <- progress::progress_bar$new(
      format = "  Simulated Annealing [:bar] :percent ETA: :eta",
      total  = niter, clear = FALSE, width = 60
    )
  }
 
  # initialize convergence check plot data.frame
  converge_tf <- 
    identical(check_converge, TRUE) || # if TRUE
    (is.numeric(check_converge) &&     # if numeric and length is 1 and is not NA and not 0
       length(check_converge) == 1 && 
       !is.na(check_converge) &&
       check_converge > 0)
  
  if (converge_tf) {
    # set vector of iterations
    check_converge <- ifelse(
      isTRUE(check_converge) || check_converge >= niter, 
      niter, check_converge
      )
    check_converge <- round(seq(1, niter, length.out = check_converge))[-1]
    
    
    non_zero_idx <- which(f_b != 0, arr.ind = TRUE)
    
    fm_iter <- 
      data.frame(
        iter    = 0, # iteration number
        phyto   = rownames(f_b)[non_zero_idx[, 1]], # phyto groups
        pigment = colnames(Fmat)[non_zero_idx[, 2]], # pigments
        ratio   = f_b[non_zero_idx] # pigment ratios
      )
  }
  
  # extract min and max F matrix thresholds and last column ratio
  wrangled <- Wrangling(f_c, min_max_mat[[1]], min_max_mat[[2]])
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
    vect <- vectorise(f_n[, -ncol(f_n)])
    oob  <- which(vect < minF | vect > maxF) 
    
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
      vect <- vectorise(f_n[, -ncol(f_n)])
      oob  <- which(vect < minF | vect > maxF) 
      
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
          "Iterations:        ", sprintf("%03d", k), "of", niter,
          "\nCurrent error:     ", round(f_c_err, 4),
          "\nNeighbour's error: ", round(f_n_err, 4),
          "\nTemperature (%):   ", round(Temp * 100, 2), "\n"
          )
        )
    }
    
    # capture f_b for convergence plot per iteration
    if (converge_tf && (k %in% check_converge)) {
      
      non_zero_idx <- which(f_b != 0, arr.ind = TRUE)
      fm_temp <- 
        data.frame(
          iter    = k, # iteration number
          phyto   = rownames(f_b)[non_zero_idx[, 1]], # phyto groups
          pigment = colnames(Fmat)[non_zero_idx[, 2]], # pigments
          ratio   = f_b[non_zero_idx] # pigment ratios
        )
      fm_iter <- rbind(fm_iter, fm_temp)
    }
  }

  final_results <- NNLS_MF_Final(f_b, S, S_Chl, S_weights)
  
  # create convergence plot
  if (converge_tf) {
    
    converge <- convergence_figure(fm_iter, niter)
    return(c(final_results, converge))
    
  }
  
  return(final_results)
}
