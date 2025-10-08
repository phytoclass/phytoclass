#' Final step for MF with prochlorococcus
#' @keywords internal
#'
#' @param Fn 
#' @param S 
#' @param S_Chl 
#' @param S_weights 
#' @param S_dvChl 
#'
#' @return
#'
#' @examples
Prochloro_NNLS_MF_Final <- function (Fn, S, S_Chl, S_weights, S_dvChl) {
  F.sum <- Prochloro_Normalise_F(Fn)[[2]]
  Fn    <- Prochloro_Normalise_F(Fn)[[1]]
  Fn    <- Fn * F.sum
  
  # determine final class abundance for non-prochloro
  Fn_wt_err <- t(Weight_error(Fn, S_weights))
  S_wt_err  <- t(Weight_error(S, S_weights))
  
  b       <- crossprod(Fn_wt_err, S_wt_err) # right hand side of linear eq
  Fn_prod <- crossprod(Fn_wt_err) # positive definite matrix with coefficients 
  
  C_new2 <- 
    RcppML::nnls(
      Fn_prod,
      b, 
      cd_maxit = 1000, 
      cd_tol = 1e-10
  )
  C_new2 <- t(C_new2)
  C_new2 <- as.matrix(C_new2)
  Cn.s2  <- rowSums(C_new2)
  Cn2    <- C_new2 / Cn.s2
  Cn2    <- as.matrix(Cn2)
  Cn2    <- Cn2 * S_Chl
  
  colnames(Cn2) <- rownames(Fn)
  colnames(Fn)  <- colnames(S)
  
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
  Cn2 <- as.data.frame(Cn2)
  
  # ---- calculate error terms ---- #
  S_residual <- S - (C_new2 %*% Fn)       # residual error
  S_rmse     <- sqrt(mean(S_residual^2))  # RMSE
  S_mae      <- colMeans(abs(S_residual)) # MAE
  
  # ---- condition number ---- #
  cd <- kappa(Fn %*% t(S))
  
  # ---- plot final results ---- #
  plt <- phyto_figure(Cn2)
  
  # row.names(Cn2) <- row.names(G)
  
  return(
    list(
      "F matrix"         = Fn, 
      "RMSE"             = S_rmse, 
      "condition number" = cd,
      "Class abundances" = Cn2, 
      "Figure"           = plt, 
      "MAE"              = S_mae, 
      "Error"            = S_residual
      )
  )
}


