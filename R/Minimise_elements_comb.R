#' Part of the steepest descent algorithm that attempts to minimize error 
#' by iteratively adjusting matrix elements
#' 
#' @keywords internal
#'
#' @param Fmat A list containing the F matrix, RMSE, and C matrix
#' @param place A vector of indices indicating which elements to adjust
#' @param S A matrix of samples (rows) and pigments (columns)
#' @param cm A vector of bounded weights for each pigment
#' @param c1_num A numeric vector (1, 2, or 3) indicating which scaler values to use
#'
#' @return A list containing the optimized F matrix with reduced error
#'
#' @examples
#' # Create sample matrices
#' F <- matrix(c(0.5, 0.3, 0.2,
#'               0.4, 0.1, 0.5), nrow=2, byrow=TRUE)
#' S <- matrix(runif(12), nrow=4)
#' cm <- c(1, 1, 1)
#' Fmat <- list(F, 0.1, matrix(1, nrow=2, ncol=4))  # F matrix, RMSE, C matrix
#' place <- c(1, 2, 3)  # elements to minimize
#' 
#' # Run optimization
#' result <- Minimise_elements_comb(Fmat, place, S, cm)
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
#' Minimise_elements <- function(Fmat, place, S, cm){   # A function that reduces every for every element that didn't reduce in index function
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
#' Minimise_elements2 <- function(Fmat, place, S, cm){   # A function that reduces every for every element that didn't reduce in index function
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
