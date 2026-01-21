#' Part of the steepest descent algorithm
#' 
#' @keywords internal
#'
#' @param Fmat The F matrix
#' @param place A vector of indices indicating which elements to adjust
#' @param S A matrix of samples (rows) and pigments (columns)
#' @param cm A vector of bounded weights for each pigment
#' @param c1_num A numeric vector (1, 2, or 3) indicating which scaler values to use
#'
#' @return A list containing the F matrix with updated ratios, the RMSE of new estimates,
#'   and an updated C matrix of estimated group contribution.
#'
#' @examples
#'  Fmat <- as.matrix(phytoclass::Fm)
#'  S <- as.matrix(phytoclass::Sm)
#'  S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
#'  place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
#'
#'  # Get F_initial from NNLS_MF as done in Steepest_Descent
#'  F_initial <- phytoclass::NNLS_MF(Fmat, S, S_weights)
#'
#'  # Run Minimise_elements_comb with c1_num = 3 (as in Steepest_Descent)
#'  result <- phytoclass:::Minimise_elements_comb(
#'    F_initial[[1]], place, S, S_weights, c1_num = 3
#'  )
Minimise_elements_comb <- function(Fmat, place, S, cm, c1_num = c(1, 2, 3)) { # A function that reduces every for every element that didn't reduce in index function
  
  f     <- Conduit(Fmat, place, S, cm, c_num = c1_num) # Calls index function
  F.new <- f[[1]] # F matrix
  n     <- f[[2]] # elements that reduce error
  if (is.null(n)) {
    n <- place
  }
  
  F.old     <- f[[3]] # old F matrix
  F.initial <- F.new # Fac_F new matrix
  
  # Fac_F new matrix
  place1 <- NULL
  if (c1_num != 1) {place1 <- place }
  
  g <- Fac_F_RR(F.new, vary = place, place = place1, S, cm, fac_rr = c1_num)
  
  if (g[[1]][[2]] < F.initial[[2]]) {
    F.new <- g[[1]]
  }
  
  n   <- g[[2]]
  res <- list(F.new, n)
  return(F.new)
  
}

# ============================================================================ #
# ---- old versions ---- #
# ============================================================================ #

#' #' Part of the steepest descent algorithm
#' #' 
#' #' @keywords internal
#' #'
#' #' @param Fmat    xx
#' #' @param place  xx
#' #' @param S   xx
#' #' @param cm  xx
#' #'
#' #' @return
#' #'
#' #' @examples
#' Minimise_elements <- function(Fmat, place, S, cm){
#'   # A function that reduces every for every element that didn't reduce
#'   # in index function
#'   f <- Conduit_3(Fmat, place, S, cm) # Calls index function
#'   F.new <- f[[1]] # F matrix
#'   n <- f[[2]] #elements that reduce error
#'   if (is.null(n)){n <- place}
#'   F.old <- f[[3]] # old F matrix
#'   F.initial <- F.new # Fac_F new matrix
#'   # Fac_F new matrix
#'   g <- Fac_F_RR3(F.new, vary = place, place, S, cm)
#'   if (g[[1]][[2]] < F.initial[[2]]){F.new <- g[[1]]}
#'   n <- g[[2]]
#'   res <- list(F.new,n)
#'   return(F.new)
#' }
#' 
#' #' A function that reduces every for every element that didn't reduce in index function
#' #'
#' #' @keywords internal
#' #'
#' #' @param Fmat   xx
#' #' @param place  xx
#' #' @param S    xx
#' #' @param cm     xx 
#' #'
#' #' @return
#' #'
#' #' @examples
#' Minimise_elements1 <- function(Fmat, place, S, cm){  
#'   # A function that reduces every for every element that didn't reduce in index function
#'   f <- Conduit_1(Fmat, place, S, cm) # Calls index function
#'   F.new <- f[[1]] # F matrix
#'   n <- f[[2]] #elements that reduce error
#'   if (is.null(n)){n <- place}
#'   F.old <- f[[3]] # old F matrix
#'   F.initial <- F.new # Fac_F new matrix
#'   # Fac_F new matrix
#'   g <-Fac_F_RR1(F.new, place, S, cm)
#'   if (g[[1]][[2]] < F.initial[[2]]){F.new <- g[[1]]}
#'   n <- g[[2]]
#'   res <- list(F.new,n)
#'   return(F.new)
#' }
#' 
#' #' A function that reduces every for every element that didn't reduce in index function
#' #' 
#' #' @keywords internal
#' #'
#' #' @param Fmat   xx
#' #' @param place  xx
#' #' @param S   xx
#' #' @param cm  xx
#' #'
#' #' @return
#' #'
#' #' @examples
#' Minimise_elements2 <- function(Fmat, place, S, cm){
#'   # A function that reduces every for every element that didn't reduce
#'   # in index function
#'   f <- Conduit_2(Fmat, place, S, cm) # Calls index function
#'   F.new <- f[[1]] # F matrix
#'   n <- f[[2]] #elements that reduce error
#'   if (is.null(n)){n <- place}
#'   F.old <- f[[3]] # old F matrix
#'   F.initial <- F.new # Fac_F new matrix
#'   # Fac_F new matrix
#'   g <-Fac_F_RR2(F.new, vary = place, place, S, cm)
#'   if (g[[1]][[2]] < F.initial[[2]]){F.new <- g[[1]]}
#'   n <- g[[2]]
#'   res <- list(F.new,n)
#'   return(F.new)
#' }
