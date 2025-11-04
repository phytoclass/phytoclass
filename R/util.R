#' Vectorise a matrix and keep non-zero elements
#' 
#' Turn each non-zero element of the F-matrix into a vector
#' 
#' @keywords internal
#'
#' @param Fmat  A matrix to vectorise
#'
#' @return A vector of non-zero pigment elements
vectorise <- function(Fmat) {
  return(Fmat[Fmat > 0])
}

#' This function normalises each column in F to row sum
#' 
#' @keywords internal
#' 
#' @param Fmat A matrix or data.frame with the last column containing non-zeros
#'
#' @return A list consisting of two components:
#'     - a matrix of pigment ratios normalized to row sums
#'     - a vector of row sums
Normalise_F <- function(Fmat) {
  Fmat  <- as.matrix(Fmat)           # convert to matrix
  F_1   <- Fmat / Fmat[, ncol(Fmat)] # divide Fmat by last column
  F.sum <- rowSums(F_1)              # sum rows
  F_1   <- F_1 / F.sum               # divide by sum
  return(list(F_1, F.sum))
}

#' Normalize F matrix specifically for Prochlorococcus pigments
#'
#' Normalizes pigment ratios differently for Prochlorococcus vs other groups,
#' using divinyl chlorophyll a for Prochlorococcus and chlorophyll a for others.
#'
#' @keywords internal
#'
#' @param Fmat Matrix with pigment ratios, where the last column is chlorophyll a
#'             and second-to-last column is divinyl chlorophyll a
#'
#' @return A list containing:
#'   \item{[[1]]}{Normalized matrix where each row is scaled by its biomass marker}
#'   \item{[[2]]}{Vector of row sums from the scaled matrix}
#'
#' @examples
#' # Create sample F matrix with Prochlorococcus
#' Fmat <- as.matrix(phytoclass::Fp)
#' result <- Prochloro_Normalise_F(Fmat)
Prochloro_Normalise_F <- function(Fmat) {
  f_new <- as.matrix(Fmat)
  n <- nrow(f_new)
  p <- ncol(f_new)
  
  # Identify Pro row (prefer rowname; else assume last row)
  pro_names <- c("pro", "prochlorococcus", "prochlorococcus-1", "pro-1")
  i_pro <- which(tolower(rownames(Fmat)) %in% pro_names)
  if (length(i_pro) != 1) i_pro <- n
  i_nonpro <- setdiff(seq_len(n), i_pro)
  
  chla   <- f_new[, p]        # last column is Chl a (Tchla)
  dvchla <- f_new[, p - 1]    # second last is dvChl a
  
  # --- Row-wise scaling so biomass pigment == 1 ---
  # Non-Pro groups: scale by Chl a of that row
  if (length(i_nonpro) > 0) {
    denom <- chla[i_nonpro]
    denom[!is.finite(denom) | denom == 0] <- 1   # guard
    f_new[i_nonpro, ] <- f_new[i_nonpro, , drop = FALSE] / denom
  }
  
  # Pro row: scale by its dvChl a
  dv_pro <- dvchla[i_pro]
  if (!is.finite(dv_pro) || dv_pro == 0) dv_pro <- 1
  f_new[i_pro, ] <- f_new[i_pro, , drop = FALSE] / dv_pro
  
  # Enforce exact biomass markers after scaling
  f_new[i_pro, p - 1] <- 1      # Pro: dvChl a == 1
  # (Chl a for Pro should already be 0; keep it as-is.)
  
  # Final step mirrors your original pipeline:
  # compute row sums AFTER pigment scaling; return the *pre-row-sum* scaled version via Fn <- Fn * F.sum
  f_sum <- rowSums(f_new)
  f_norm <- f_new / f_sum
  return(list(as.matrix(f_norm), f_sum))
}

#' Normalise matrix to row sum
#' 
#' This function normalises each column in S to row sum
#' 
#' @keywords internal
#'
#' @param S A matrix or data.frame to be normalized
#'
#' @return A matrix
#'
#' @examples
#' # Create a sample matrix
#' S <- as.matrix(phytoclass::Sm)
#' normalized <- phytoclass:::Normalise_S(S)
#' 
Normalise_S <- function(S){
  # Normalise to unit row sum
  S <- as.matrix(S)
  S <- S / rowSums(S)
  return(S)
}

#' Add weights to the data, bound at a maximum.
#' 
#' @param S   Sample data matrix – a matrix of pigment samples
#' @param weight.upper.bound  Upper bound for weights (default is 30)        
#'
#' @return A vector with upper bounds for weights
#' @export
#'
#' @examples
#' Bounded_weights(Sm, weight.upper.bound = 30)
#' 
Bounded_weights <- function(S, weight.upper.bound = 30) {
  n <- colMeans(S)
  S <- n^-1
  S[S > weight.upper.bound] <- weight.upper.bound
  S[length(S)] <- 1
  return(S)
}

#' Wrangle data to vectors
#' 
#' Converts data-types and selects data for randomisation in 
#' the simulated annealing algorithm
#' 
#' @keywords internal
#'
#' @param Fl      A matrix of the initial F matrix (i.e. pigment ratio matrix)
#' @param min.val A vector of the minimum values for each non-zero pigment ratios 
#' @param max.val A vector of the maximum values for each non-zero pigment ratios
#'
#' @return
#'     A list containing following components:
#'     - A vector Fmin with the minimum pigment ratio values
#'     - A vector Fmax with the maximum pigment ratio values
#'     - A vector SE with the current pigment ratio values
#'     - A vector chlv with the pigment ratio values for the last column in Fl
Wrangling <- function(Fl, min.val, max.val) {
  
  # set up initial F, Fmin, Fmax, matrix by removing Tchla column
  Fd <- Fmin <- Fmax <- as.matrix(Fl)[, -ncol(Fl)]
  
  # set all non-zero elements of F matrix to the minimum and maximum values
  Fmin[Fmin > 0] <- min.val
  Fmax[Fmax > 0] <- max.val
  
  # extract chlorophyll-a values once weighted to rowsums for initial F matrix
  chlv <- Fl[, ncol(Fl)]
  
  # multiply the minimum value by weighted chlorophyll to updated ratios
  Fmin <- Fmin * chlv
  Fmax <- Fmax * chlv
  
  # vectorise function outputs all non-zero elements as a vector (excluding chl column)
  Fmin <- vectorise(Fmin)
  Fmax <- vectorise(Fmax)
  SE   <- vectorise(Fd)
  
  return(list(Fmin, Fmax, SE, chlv))
}


#' Prochloro Wrangling
#' 
#' @keywords internal
#'
#' @param Fl Initial F matrix (i.e. pigment ratio matrix)
#' @param min.val A vector of minimum values for each non-zero pigment ratio
#' @param max.val A vector of maximum values for each non-zero pigment ratio
#'
#' @return A list containing vectorized min/max bounds and current values for
#'   Prochlorococcus-aware normalization
Prochloro_Wrangling <- function(Fl, min.val, max.val) {
  
  # set up initial F, Fmin, Fmax, matrix by removing Tchla column
  Fd <- Fl
  Fmin <- Fmax <- as.matrix(Fl)[, -ncol(Fl)]
  
  Fmin[Fmin > 0] <- min.val
  Fmax[Fmax > 0] <- max.val
  
  chlv  <- Fl[, ncol(Fl)]      # Chl a
  chlvp <- Fl[, ncol(Fl) - 1]  # dvChl a
  chlep <- Fl[, ncol(Fl) - 1]
  chlep[length(chlep)] <- 1
  
  # Non-Prochloro rows scaled by Chl a (exclude dvChl column)
  Fmin[-nrow(Fmin), -ncol(Fmin)] <- Fmin[-nrow(Fmin), -ncol(Fmin)] * chlv[-length(chlv)]
  Fmax[-nrow(Fmax), -ncol(Fmax)] <- Fmax[-nrow(Fmax), -ncol(Fmax)] * chlv[-length(chlv)]
  
  # Prochlorococcus row scaled by dvChl a
  Fmin[nrow(Fmin), ] <- Fmin[nrow(Fmin), ] * chlvp[length(chlvp)]
  Fmax[nrow(Fmax), ] <- Fmax[nrow(Fmax), ] * chlvp[length(chlvp)]
  
  # Vectorised CORE (base-R column-major) and remove zeros
  Fmin_vec <- vectorise(Fmin)
  Fmax_vec <- vectorise(Fmax)
  SE_vec   <- vectorise(Fl[, 1:(ncol(Fmin) - 2)])
  
  return(list(Fmin_vec, Fmax_vec, SE_vec, chlv, chlvp))
}

#' Apply weights to F/S matrices by diagonal multiplication
#' 
#' @keywords internal
#'
#' @param S Matrix to be weighted
#' @param cm Vector of weights to be applied to columns of S
#'
#' @return A matrix with weighted columns (S %*% diag(cm))
#'
#' @examples
#' # Create sample matrix and weights
#'  Fmat <- as.matrix(phytoclass::Fm)
#'  S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
#' weighted <- Weight_error(Fmat, S_weights)
Weight_error <- function(S, cm){
  S <- S %*% diag(cm)
  return(S)
}

#' Calculate the mean condition number for randomized F matrices
#' 
#' Performs multiple simulations with randomized F matrices within given bounds
#' to assess the numerical stability of the system.
#' 
#' @keywords internal
#'
#' @param S Sample matrix of pigment measurements
#' @param Fn Initial F matrix of pigment ratios
#' @param min.val Optional vector of minimum values for each non-zero pigment ratio
#' @param max.val Optional vector of maximum values for each non-zero pigment ratio
#'
#' @return Numeric value representing the mean condition number from 1000 simulations
#'
#' @examples
#' # Create sample matrices
#' Fmat <- as.matrix(phytoclass::Fm)
#' S <- as.matrix(phytoclass::Sm)
#' min_max <- phytoclass::min_max
#' min.val <- min_max[[3]]
#' max.val <- min_max[[4]]
#' # Use only non-chla columns for Fmat and S, as in simulated_annealing
#' Fmat_sub <- Fmat[, -ncol(Fmat)]
#' S_sub <- S[, -ncol(S)]
#' # Calculate mean condition number
#' cond <- phytoclass:::Condition_test(S_sub, Fmat_sub, min.val, max.val)
Condition_test <- function(S, Fn, min.val = NULL, max.val = NULL) {
  if (is.null(min.val) & is.null(max.val)) {
    min_max <- Default_min_max(phytoclass::min_max, Fn)
    min.val <- min_max[[1]]
    max.val <- min_max[[2]]
  }
  condition_number <- function(S, f_mat, min.val, max.val) {
    f_non_zero <- length(vectorise(f_mat))
    rand       <- vector(length = f_non_zero)
    for (i in seq(f_non_zero)) {
      rand[i] <- stats::runif(1, min = min.val[i], max = max.val[i])
    }
    f_mat[f_mat > 0] <- rand
    return(kappa(f_mat %*% t(S)))
  }
  sn <- replicate(n = 1000, condition_number(S, Fn, min.val, max.val))
  return(mean(sn))
}

