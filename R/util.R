#' Vectorise a matrix and keep non-zero elements
#' 
#' Turn each non-zero element of the F-matrix into a vector
#' 
#' @keywords internal
#'
#' @param Fmat  A matrix to vectorise
#'
#' @return A vector of non-zero pigment elements
#'
#' @examples
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
#'
#' @examples
Normalise_F <- function(Fmat) {
  Fmat  <- as.matrix(Fmat)           # convert to matrix
  F_1   <- Fmat / Fmat[, ncol(Fmat)] # divide Fmat by last column
  F.sum <- rowSums(F_1)              # sum rows
  F_1   <- F_1 / F.sum               # divide by sum
  return(list(F_1, F.sum))
}

#' Normalise F for prochloro
#'
#' @keywords internal
#'
#' @param Fmat 
#'
#' @return
#'
#' @examples
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
#'
#' @keywords internal
#'
#' @param S   A matrix or data.frame
#'
#' @return A matrix
#'
#' @examples
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
#' @examples
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
#' @param Fl 
#' @param min.val 
#' @param max.val 
#'
#' @return
#'
#' @examples
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

#' Apply weights to F/S matrices
#' 
#' @keywords internal
#'
#' @param S  xx
#' @param cm  xx
#'
#' @return A matrix
#'
#' @examples
Weight_error <- function(S, cm){
  S <- S %*% diag(cm)
  return(S)
}