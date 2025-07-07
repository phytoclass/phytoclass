#' Performs the non-negative matrix factorisation for given phytoplankton 
#' pigments and pigment ratios, to attain an estimate of phytoplankton 
#' class abundances.
#'
#' @param Fn   Pigment to Chl a matrix
#' @param S   Sample data matrix – a matrix of pigment samples
#' @param S_weights  Weights for each column
#'
#' @return A list containing 
#' \enumerate{
#'  \item The F matrix (pigment: Chl *a*) ratios
#'  \item The root mean square error (RMSE)
#'  \item The C matrix (class abundances for each group)
#'  }
#' @export
#'
#' @examples
#' MC <- Matrix_checks(Sm,Fm)
#' Snew <- MC$Snew
#' Fnew <- MC$Fnew
#' S_weights <- Bounded_weights(Snew, weight.upper.bound = 30)
#' NNLS_MF(Fnew, Snew, S_weights)
#'
NNLS_MF <- function(Fn, S, S_weights = NULL) {
  if (is.null(S_weights)) {
    S_weights <- as.vector(rep(1, ncol(S)))
  }

  Fn_wt_err <- t(Weight_error(Fn, S_weights))
  S_wt_err  <- t(Weight_error(S, S_weights))
  
  b       <- crossprod(Fn_wt_err, S_wt_err)
  Fn_prod <- crossprod(Fn_wt_err)
  
  C_new2  <- 
    RcppML::nnls(
      Fn_prod,
      b,
      cd_maxit = 1000, 
      cd_tol   = 1e-8
    )
  C_new2 <- t(C_new2)
  Cn2    <- Normalise_S(C_new2) # Row sums to one
  
  # b <- crossprod(t(Weight_error(Fn, S_weights)), t(Weight_error(S, S_weights)))
  # C_new2 <- t(RcppML::nnls(crossprod(t(Weight_error(Fn, S_weights))),
  #   b,
  #   cd_maxit = 1000, cd_tol = 1e-8
  # ))
  # C_new2 <- as.matrix(C_new2)
  # Cn.s2  <- rowSums(C_new2)
  # Cn2    <- C_new2 / Cn.s2 # Row sums to one
  # Cn2    <- as.matrix(Cn2)
  colnames(Cn2) <- rownames(Fn)
  
  # error <- Metrics::rmse(S,(C_new2%*%Fn))
  error <- sqrt(mean((S - C_new2 %*% Fn)^2)) # RMSE
  
  return(list("F matrix " = Fn, "RMSE" = error, "C matrix" = Cn2))
}



