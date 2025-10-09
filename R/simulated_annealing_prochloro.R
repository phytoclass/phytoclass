#' Perform simulated annealing algorithm for samples with divinyl chlorophyll and prochlorococcus. 
#' Chlorophyll must be the final column of both S and F matrices, with Divinyl Chlorophyll a the 2nd to last column. 
#' See how the example Sp and Fp matrices are organised.
#'
#' @param S   Sample data matrix – a matrix of pigment samples
#' @param Fmat   Pigment to Chl a matrix
#' @param user_defined_min_max data frame with some format as min_max built-in data
#' @param do_matrix_checks     This should only be set to TRUE when using the default values. This will remove pigment columns that have column sums of 0. Set to FALSE if using customised names for pigments and phytoplankton groups
#' @param niter Number of iterations (default is 500)
#' @param step  Step ratio used (default is 0.009)
#' @param weight.upper.bound Upper limit of the weights applied (default value is 30). 
#' @param verbose Logical value. Output error and temperature at each iteration. Default value of TRUE
#' @param seed Set seed number to reproduce the same results
#' @param check_converge  TRUE/FALSE/integer; set the number of F matrices to 
#'                        for convergence checking
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
simulated_annealing_Prochloro <- function(
    S, 
    Fmat                 = NULL,
    user_defined_min_max = NULL,
    do_matrix_checks     = TRUE, 
    niter                = 500, 
    step                 = 0.009, 
    weight.upper.bound   = 30, 
    verbose              = TRUE,
    seed                 = NULL,
    check_converge       = 100
  ) {
  
  # set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (is.null(Fmat)) {
    Fmat <- phytoclass::Fp
  }
  
  # Strip text cols if any
  if (is.data.frame(S)) {
    char_cols <- sapply(S, is.character)
    S         <- S[, !char_cols]
  }
  
  if (!is.matrix(S)) {
    S <- as.matrix(S)
  }
  
  # Optional matrix checks / alignment
  if (do_matrix_checks) {
    mat_n <- Matrix_checks(as.matrix(S), as.matrix(Fmat))
    S     <- as.matrix(mat_n[[1]])
    Fmat  <- as.matrix(mat_n[[2]])
  }
  
  # Save Chl & dvChl before normalising S (used at the very end)
  S_Chl   <- S[, ncol(S)]  # chla 
  S_dvChl <- S[, ncol(S) - 1]  # dvchla 
  
  # Normalise S and build weights
  S  <- Normalise_S(S)
  S_weights <- Bounded_weights(S, weight.upper.bound)
  
  # SAALS compatibility: 'place' includes dvChl, excludes Chl (as in your original)
  place_full     <- which(Fmat[, -ncol(Fmat)] > 0) # non-zero, non-chla pigments
  
  # min/max lookup for SAALS domain
  if (is.null(user_defined_min_max)) {
    min_max_mat <- Default_min_max(phytoclass::min_max, Fmat[, -ncol(Fmat)], place_full)
  } else {
    min_max_mat <- Default_min_max(user_defined_min_max, Fmat[, -ncol(Fmat)], place_full)
  }
  
  # Condition number (unchanged logic)
  condition.test <- Condition_test(
    S[,    -ncol(S)],
    Fmat[, -ncol(Fmat)],
    # min.val, max.val
    min_max_mat[[1]], min_max_mat[[2]]
  )
  if (verbose) {
    message(paste0("\nCondition number = ", round(condition.test), "\n\n"))
  }
  if (condition.test > 1e5) {
    print("Abort process: condition number of S matrix greater than 100 000\n")
  }
  
  # ---- start iteration process ---- #
  
  # Fixed binary mask and initial NNLS
  Fmat <- ifelse(Fmat > 0, 1, 0)
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
  
  # ---- Static bounds for the CORE (exclude dvChl & Chl), computed once from f_c
  W0        <- Prochloro_Wrangling(f_c, min_max_mat[[1]], min_max_mat[[2]])  # returns vectorised core min/max
  minF_fix  <- W0[[1]]
  maxF_fix  <- W0[[2]]
  # ---------------------------------------------------------------------------
  
  # # Helpers
  # pick_best <- function(D, fallback) {
  #   if (length(D) == 0) return(fallback)
  #   Dn <- suppressWarnings(vapply(D, function(i) i[[2]], numeric(1)))
  #   ok <- is.finite(Dn)
  #   if (!any(ok)) return(fallback)
  #   D[[ which(ok)[ which.min(Dn[ok]) ] ]]
  # }
  # safe_bounds <- function(a,b) {
  #   all(is.finite(a)) && all(is.finite(b)) && all(a <= b, na.rm = TRUE)
  # }
  
  step <- 1 - step
  
  for (k in 1:niter) {
    
    if (!verbose) pb$tick()
    
    Temp    <- step^k
    
    # Always source dvChl/Chl from the same matrix we perturb
    chlv  <- f_c[, ncol(f_c)]        # Chl a (Tchla)
    chlvp <- f_c[, ncol(f_c) - 1]    # dvChl a
    
    # Seed neighbour (all-at-once)
    seed_nb <- Prochloro_Random_Neighbour_2(
      f_c, Temp, chlv, f_c, place_full, S, S_weights,
      minF_fix, maxF_fix, chlvp, Fi_mask = Fmat
    )
    
    # Candidate pool
    num_loop <- ifelse(k > niter - 20, 300, 120)
    
    # set list to store random neighbors
    rand_itr_err <- rand_itr <- vector("list", num_loop) 
    
    for (i in seq(num_loop)) {
      chlv  <- f_c[, ncol(f_c)]
      chlvp <- f_c[, ncol(f_c) - 1]
      temp_rand <- Prochloro_Random_Neighbour_2(
        f_c, Temp, chlv, f_c, place_full, S, S_weights, 
        minF_fix, maxF_fix, chlvp, 
        Fi_mask = Fmat
      )
      
      rand_itr[[i]]     <- temp_rand
      rand_itr_err[[i]] <- temp_rand[[2]] # extract RMSE
      
    }

    # select neighbor with lowest RMSE
    # TODO: add seed as fallback
    low_indx      <- which.min(rand_itr_err)
    new_neighbour <- rand_itr[[low_indx]]
    
    # steepest descent
    num_loop2     <- ifelse(Temp > 0.3, 10, 2)
    new_neighbour <- Steepest_Descent(new_neighbour[[1]], place_full, S, S_weights, 
                                      num.loops = num_loop2)
    
    f_n     <- new_neighbour[[1]]
    
    # ========= delete start
    # nrep <- if (k > niter - 20) 300 else 120
    # D <- vector("list", nrep)
    # for (i in seq_len(nrep)) {
    #   chlv  <- f_c[, ncol(f_c)]
    #   chlvp <- f_c[, ncol(f_c) - 1]
    #   D[[i]] <- Prochloro_Random_Neighbour_2(
    #     f_c, Temp, chlv, f_c, place_full, S, S_weights,
    #     minF_fix, maxF_fix, chlvp, Fi_mask = Fmat
    #   )
    # }
    # new_neighbour <- pick_best(D, seed_nb)
    # 
    # # all(is.finite(a)) && all(is.finite(b)) && all(a <= b, na.rm = TRUE)
    # # Local search (unchanged API)
    # if (safe_bounds(min_max_mat[[1]], min_max_mat[[2]])) {
    #   if (Temp > 0.3) {
    #     new_neighbour <- SAALS(new_neighbour[[1]], min_max_mat[[1]], min_max_mat[[2]], place_full, S, S_weights, num.loops = 10)
    #   } else {
    #     new_neighbour <- SAALS(new_neighbour[[1]], min_max_mat[[1]], min_max_mat[[2]], place_full, S, S_weights, num.loops = 2)
    #   }
    # }
    # ========= delete end
    
    

    # f_n_err <- new_neighbour[[2]]
    
    # # Enforce static CORE bounds by projection if needed (dvChl & Chl excluded)
    # core_vec_n <- as.vector(f_n[, 1:(ncol(f_n) - 2)])
    # bad        <- which(core_vec_n < minF_fix | core_vec_n > maxF_fix)
    # if (length(bad) > 0) {
    #   core_vec_n[bad] <- pmin(pmax(core_vec_n[bad], minF_fix[bad]), maxF_fix[bad])
    #   core_n          <- f_n[, 1:(ncol(f_n) - 2)]
    #   core_n[]        <- core_vec_n
    #   f_n             <- cbind(core_n, f_n[, ncol(f_n) - 1], f_n[, ncol(f_n)])
    #   colnames(f_n)   <- colnames(S)
    #   tmp             <- NNLS_MF(f_n, S, S_weights)
    #   f_n             <- tmp[[1]]
    #   f_n_err         <- tmp[[2]]
    # }
    
    # check if ratios are out of bounds (min\max)
    vect    <- vectorise(f_n[, seq(ncol(f_n) - 2)])
    oob     <- which(vect < min_max_mat[[1]][-length(min_max_mat[[1]])]
                     | vect >  min_max_mat[[2]][-length( min_max_mat[[2]])]
                     )
    
    while (length(oob) > 0) {
      
      oob_indx  <- place_full[oob] # where in F matrix is the ratio out of bounds
      num_loop  <- ifelse(k > niter - 20, 300, 120)
      
      # set list to store random neighbors
      rand_itr2_err <- rand_itr2 <- vector("list", num_loop)
      
      for (i in seq(num_loop)) {
        chlv  <- f_c[, ncol(f_c)]
        chlvp <- f_c[, ncol(f_c) - 1]
        
        temp_rand <- Prochloro_Random_Neighbour(
          f_n, Temp, chlv, chlvp, oob_indx, S, S_weights, 
          minF_fix, maxF_fix
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
      oob     <- which(vect < min_max_mat[[1]][-length(min_max_mat[[1]])]
                       | vect >  min_max_mat[[2]][-length( min_max_mat[[2]])]
                       )
    }
    
    
    # -------------------- ORIGINAL ACCEPTANCE RULE (wrapped in finite guard) --------------------
    # check RMSE of neighbor is better than current 
    f_n_err  <- new_neighbour[[2]] # RMSE
    if (is.finite(f_n_err) && is.finite(f_c_err) &&
        (f_n_err < f_c_err || 
         exp(-(f_n_err - f_c_err)) < stats::runif(1, 0, 1))
        ) {
      f_c     <- f_n
      f_c_err <- f_n_err
    }
    # -------------------------------------------------------------------------------------------
    
    # update F matrix and RMSE when neighbor F matrix has lower RMSE
    if (is.finite(f_n_err) && f_n_err < f_b_err) {
      f_b <- f_n
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
  
  # Final NNLS + aggregation (unchanged)
  final_results <- Prochloro_NNLS_MF_Final(f_b, S, S_Chl, S_weights, S_dvChl)
  
  # create convergence plot
  if (converge_tf) {

    converge <- convergence_figure(fm_iter, niter)
    return(c(final_results, converge))

  }
  
  return(final_results)
}

