#' Perform matrix factorisation for phytoplankton pigments and pigments ratios
#' 
#' Performs the non-negative matrix factorisation for given phytoplankton 
#' pigments and pigment ratios, to attain an estimate of phytoplankton 
#' class abundances.
#' 
#' Unlike NNLS_ML(), it also removes any weighting and normalisation, and 
#' also multiplies relative abundances by chlorophyll values to determine
#' the biomass of phytoplankton groups.
#' 
#' @keywords internal
#'
#' @param Fn xx
#' @param S   xx
#' @param S_Chl   xx
#' @param S_weights  xx
#' @param S_dvChl xx
#'
#' @return
#'
#' @examples  
NNLS_MF_Final <- function(Fn, S, S_Chl, S_weights, S_dvChl = NULL) {
  check_pro <- any(tolower(colnames(Fn)) %in% c("dvchl", "dvchla", "chlvp"))
  if (check_pro) {
    F_norm <- Prochloro_Normalise_F(Fn)
  } else {
    F_norm <- Normalise_F(Fn)
  }
  # normalize F matrix
  # F_norm <- Normalise_F(Fn)
  Fn     <- F_norm[[1]] * F_norm[[2]]

  Fn_wt_err <- t(Weight_error(Fn, S_weights))
  S_wt_err  <- t(Weight_error(S, S_weights))
  
  b       <- crossprod(Fn_wt_err, S_wt_err) # right hand side of linear eq
  Fn_prod <- crossprod(Fn_wt_err) # positive definite matrix with coefficients 

  # ---- calc NNLS ---- #
  C_new2  <- 
    RcppML::nnls(
      Fn_prod,
      b,
      cd_maxit = 1000, 
      cd_tol   = 1e-10
    )
  C_new2 <- t(C_new2)

  C_new2        <- as.matrix(C_new2)
  Cn.s2         <- rowSums(C_new2)
  Cn.s2         <- ifelse(Cn.s2 == 0, 1, Cn.s2)
  Cn2           <- C_new2 / Cn.s2
  Cn2           <- as.matrix(Cn2)
  Cn2           <- Cn2 * S_Chl
  colnames(Cn2) <- rownames(Fn)
  colnames(Fn)  <- colnames(S)
  Cn2           <- as.data.frame(Cn2)
  rownames(Cn2) <- rownames(S)
  
  if (check_pro) {
    # determine final class abundance for prochloro
    Fn_wt_err2 <- t(Weight_error(Fn[nrow(Fn), -ncol(Fn)], S_weights[-length(S_weights)]))
    S_wt_err2  <- t(Weight_error(S[, -ncol(S)], S_weights[-length(S_weights)]))

    Pb       <- crossprod(Fn_wt_err2, S_wt_err2)
    Fn_prod2 <- crossprod(Fn_wt_err2)

    PC_new2 <-
      RcppML::nnls(
        Fn_prod2,
        Pb,
        cd_maxit = 1000,
        cd_tol   = 1e-10
      )
    PC_new2 <- t(PC_new2)
    PCn.s2  <- rowSums(PC_new2)
    PCn2    <- PC_new2 / PCn.s2
    PCn2    <- PCn2 * S_dvChl

    Cn2[, ncol(Cn2)] <- PCn2
  }
  
  # ---- calculate error terms ---- #
  S_residual <- S - (C_new2 %*% Fn)       # residual error
  S_rmse     <- sqrt(mean(S_residual^2))  # RMSE
  S_mae      <- colMeans(abs(S_residual)) # MAE
  
  # ---- condition number ---- #
  cd <- kappa(Fn %*% t(S))
  
  # ---- plot final results ---- #
  plt <- phyto_figure(Cn2)

  return(list(
    "F matrix"         = Fn,
    "RMSE"             = S_rmse,
    "condition number" = cd,
    "Class abundances" = Cn2,
    "Figure"           = plt,
    "MAE"              = S_mae,
    "Error"            = S_residual
  )) 
}
