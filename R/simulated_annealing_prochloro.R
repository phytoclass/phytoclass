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
simulated_annealing_Prochloro <- function (
    S, 
    Fmat = NULL,
    user_defined_min_max = NULL,
    do_matrix_checks = TRUE, 
    niter = 500, 
    step = 0.009, 
    weight.upper.bound = 30, 
    verbose = TRUE
) {
  if (is.null(Fmat)) Fmat <- phytoclass::Fp
  
  # Strip text cols if any
  if (is.data.frame(S)) {
    char_cols <- sapply(S, is.character)
    S <- S[, !char_cols]
  }
  
  # Optional matrix checks / alignment
  if (do_matrix_checks) {
    L <- Matrix_checks(as.matrix(S), as.matrix(Fmat))
    S    <- as.matrix(L[[1]])
    Fmat <- as.matrix(L[[2]])
  }
  
  # Save Chl & dvChl before normalising S (used at the very end)
  S_Chl   <- S[, ncol(S)]
  S_dvChl <- S[, ncol(S) - 1]
  
  # Normalise S and build weights
  S  <- Normalise_S(S)
  cm <- Bounded_weights(S, weight.upper.bound)
  
  # SAALS compatibility: 'place' includes dvChl, excludes Chl (as in your original)
  place_full <- which(Fmat[, 1:(ncol(Fmat) - 1)] > 0)
  
  # min/max lookup for SAALS domain
  if (is.null(user_defined_min_max)) {
    K <- Default_min_max(phytoclass::min_max, Fmat[, 1:(ncol(Fmat) - 1)], place_full)
  } else {
    K <- Default_min_max(user_defined_min_max,     Fmat[, 1:(ncol(Fmat) - 1)], place_full)
  }
  min.val <- K[[1]]
  max.val <- K[[2]]
  
  # Condition number (unchanged logic)
  condition.test <- Condition_test(
    S[,    1:(ncol(S)    - 1)],
    Fmat[, 1:(ncol(Fmat) - 1)],
    min.val, max.val
  )
  if (verbose) {
    message(paste0("\nCondition number = ", round(condition.test), "\n\n"))
  }
  if (condition.test > 1e5) {
    print("Abort process: condition number of S matrix greater than 100 000\n")
  }
  
  # Fixed binary mask and initial NNLS
  Fi <- ifelse(Fmat > 0, 1, 0)
  nc <- NNLS_MF(Fi, S, cm)
  s_b <- s_c <- s_n <- nc[[1]]
  f_b <- f_c <- f_n <- nc[[2]]
  
  # initialize progress bar if verbose is FALSE
  if (!verbose) {
    pb <- progress::progress_bar$new(
      format = "  Simulated Annealing [:bar] :percent ETA: :eta",
      total  = niter, clear = FALSE, width = 60
    )
  }
  
  
  # ---- Static bounds for the CORE (exclude dvChl & Chl), computed once from s_c
  W0        <- Prochloro_Wrangling(s_c, min.val, max.val)  # returns vectorised core min/max
  minF_fix  <- W0[[1]]
  maxF_fix  <- W0[[2]]
  # ---------------------------------------------------------------------------
  
  # Helpers
  pick_best <- function(D, fallback) {
    if (length(D) == 0) return(fallback)
    Dn <- suppressWarnings(vapply(D, function(i) i[[2]], numeric(1)))
    ok <- is.finite(Dn)
    if (!any(ok)) return(fallback)
    D[[ which(ok)[ which.min(Dn[ok]) ] ]]
  }
  safe_bounds <- function(a,b) {
    all(is.finite(a)) && all(is.finite(b)) && all(a <= b, na.rm = TRUE)
  }
  
  for (k in 1:niter) {
    
    if (!verbose) pb$tick()
    
    Temp <- (1 - step)^k
    
    # Always source dvChl/Chl from the same matrix we perturb
    chlv  <- s_c[, ncol(s_c)]        # Chl a (Tchla)
    chlvp <- s_c[, ncol(s_c) - 1]    # dvChl a
    
    # Seed neighbour (all-at-once)
    seed_nb <- Prochloro_Random_Neighbour_2(
      s_c, Temp, chlv, s_c, place_full, S, cm,
      minF_fix, maxF_fix, chlvp, Fi_mask = Fi
    )
    
    # Candidate pool
    nrep <- if (k > niter - 20) 300 else 120
    D <- vector("list", nrep)
    for (i in seq_len(nrep)) {
      chlv  <- s_c[, ncol(s_c)]
      chlvp <- s_c[, ncol(s_c) - 1]
      D[[i]] <- Prochloro_Random_Neighbour_2(
        s_c, Temp, chlv, s_c, place_full, S, cm,
        minF_fix, maxF_fix, chlvp, Fi_mask = Fi
      )
    }
    new_neighbour <- pick_best(D, seed_nb)
    
    # Local search (unchanged API)
    if (safe_bounds(min.val, max.val)) {
      if (Temp > 0.3) {
        new_neighbour <- SAALS(new_neighbour[[1]], min.val, max.val, place_full, S, cm, num.loops = 10)
      } else {
        new_neighbour <- SAALS(new_neighbour[[1]], min.val, max.val, place_full, S, cm, num.loops = 2)
      }
    }
    
    s_n <- new_neighbour[[1]]
    f_n <- new_neighbour[[2]]
    
    # Enforce static CORE bounds by projection if needed (dvChl & Chl excluded)
    core_vec_n <- as.vector(s_n[, 1:(ncol(s_n) - 2)])
    bad <- which(core_vec_n < minF_fix | core_vec_n > maxF_fix)
    if (length(bad) > 0) {
      core_vec_n[bad] <- pmin(pmax(core_vec_n[bad], minF_fix[bad]), maxF_fix[bad])
      core_n <- s_n[, 1:(ncol(s_n) - 2)]
      core_n[] <- core_vec_n
      s_n <- cbind(core_n, s_n[, ncol(s_n) - 1], s_n[, ncol(s_n)])
      colnames(s_n) <- colnames(S)
      tmp <- NNLS_MF(s_n, S, cm)
      s_n <- tmp[[1]]
      f_n <- tmp[[2]]
    }
    
    # -------------------- ORIGINAL ACCEPTANCE RULE (wrapped in finite guard) --------------------
    if (is.finite(f_n) && is.finite(f_c) &&
        (f_n < f_c || exp(-(f_n - f_c)) < stats::runif(1, 0, 1))) {
      s_c <- s_n
      f_c <- f_n
    }
    # -------------------------------------------------------------------------------------------
    
    if (verbose) {
      message(paste("Current error: ", round(f_c, 4)))
      message(paste("Neighbour's error: ", round(f_n, 4)))
      message(paste("Temperature (%): ", round(Temp * 100, 2)))
      message(" ")
    }
    
    if (is.finite(f_n) && f_n < f_b) {
      s_b <- s_n
      f_b <- f_n
    }
  }
  
  # Final NNLS + aggregation (unchanged)
  return(Prochloro_NNLS_MF_Final(s_b, S, S_Chl, cm, S_dvChl))
}

