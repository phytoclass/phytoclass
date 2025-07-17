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
#'
#' @return
#'
#' @examples  
NNLS_MF_Final <- function(Fn, S, S_Chl, S_weights){
  # normalize F matrix
  F_norm <- Normalise_F(Fn)
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
  Cn.s2_safe    <- ifelse(Cn.s2 == 0, 1, Cn.s2)
  Cn2           <- C_new2 / Cn.s2_safe
  Cn2           <- as.matrix(Cn2)
  Cn2           <- Cn2 * S_Chl
  colnames(Cn2) <- rownames(Fn)
  colnames(Fn)  <- colnames(S)
  Cn2           <- as.data.frame(Cn2)
  rownames(Cn2) <- rownames(S)
  
  # ---- calculate error terms ---- #
  S_rmse     <- sqrt(mean((S - C_new2 %*% Fn)^2))  # RMSE
  S_residual <- S - (C_new2 %*% Fn)                # residual error
  S_mae      <- colMeans(abs((C_new2 %*% Fn) - S)) # MAE
  
  # ---- condition number ---- #
  cd <- kappa(Fn %*% t(S))
  
  # NULL assignment to stop NOTE during the package "Check"
  #  -  no visible binding for global variable
  vals    <- NULL
  row_num <- NULL
  
  # ---- plot final results ---- #
  colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#009E73",
                         "#001E73", "#013E73")

  PLE  <- tidyr::pivot_longer(
    data      = cbind(Cn2, "row_num" = seq(nrow(Cn2))), 
    cols      = -row_num, 
    names_to  = 'names', 
    values_to = 'vals'
    )
  
  plt <- 
    ggplot2::ggplot(PLE, ggplot2::aes_string(x = row_num, y = vals, fill = names)) +
    ggplot2::geom_area() +
    ggplot2::scale_color_manual(values = colorBlindGrey8) +
    ggplot2::scale_fill_manual(values = colorBlindGrey8) +
    ggplot2::xlab("Sample number") +
    ggplot2::ylab("Chl a concentrations") +
    ggplot2::theme_bw()

  
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
