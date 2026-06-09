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
#'  Fmat <- as.matrix(phytoclass::Fm)
#'  S <- as.matrix(phytoclass::Sm)
#'  S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
#'  place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
#'  num.loops <- 2
#'  # Run Steepest_Descent
#'  result <- phytoclass:::Steepest_Descent(Fmat, place, S, S_weights, num.loops)
#'
NNLS_MF <- function(Fn, S, S_weights = NULL) {
  if (is.null(S_weights)) {
    S_weights <- as.vector(rep(1, ncol(S)))
  }
  
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
      cd_tol   = 1e-8
    )
  # TODO: add in when RcppML updates in CRAN to v1.0.0, then require this 
  # version of RcppML
  # C_new2  <- RcppML::nnls(w = Fn_prod, A = b, cd_maxit = 1000, cd_tol = 1e-8)
  
  C_new2        <- t(C_new2)
  Cn2           <- Normalise_S(C_new2) # Row sums to one
  colnames(Cn2) <- rownames(Fn)
  
  # ---- calculate error term ---- #
  error <- sqrt(mean((S - C_new2 %*% Fn)^2)) # RMSE
  
  return(list("F matrix " = Fn, "RMSE" = error, "C matrix" = Cn2))
}

#' Per-sample equality + non-negativity constrained least squares via Lawson-Hanson
#'
#' @param Fn k × p  (class × pigment F matrix)
#' @param S  n × p  (sample × pigment, optionally normalised)
#' @param S_weights  pigment weights for each column
#' @param equality_sum: NULL, a scalar, or a length-n vector for per-row 
#'         equality target. If NULL, only the non-negativity constraint applies.
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
#'  Fmat <- as.matrix(phytoclass::Fm)
#'  S <- as.matrix(phytoclass::Sm)
#'  S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
#'  place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
#'  num.loops <- 2
#'  # Run Steepest_Descent
#'  result <- phytoclass:::Steepest_Descent(Fmat, place, S, S_weights, num.loops)
#'
nnls_lsei <- function(Fn, S, S_weights, equality_sum = NULL) {
  
  Fn_w <- t(Weight_error(Fn, S_weights))
  S_w  <- Weight_error(S, S_weights)
  k    <- nrow(Fn)
  n    <- nrow(S)
  C    <- matrix(0, nrow = n, ncol = k,
                 dimnames = list(rownames(S), rownames(Fn)))
  
  if (!is.null(equality_sum) && length(equality_sum) == 1L) {
    equality_sum <- rep(equality_sum, n)
  }
  
  for (i in seq_len(n)) {
    sol <- tryCatch(
      if (is.null(equality_sum)) {
        limSolve::lsei(
          A = Fn_w, B = S_w[i, ],
          G = diag(k), H = rep(0, k),
          verbose = FALSE)
      } else {
        limSolve::lsei(
          A = Fn_w, B = S_w[i, ],
          E = matrix(1, nrow = 1, ncol = k), 
          F = equality_sum[i],
          G = diag(k), H = rep(0, k),
          verbose = FALSE)
      },
      error = function(e) list(X = rep(0, k))
    )
    C[i, ] <- sol$X
  }
  C
}

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
#' @param Fn F matrix with pigment ratios for each phytoplankton class
#' @param S Sample data matrix of pigment measurements
#' @param S_Chl Vector of chlorophyll a concentrations for each sample
#' @param S_weights Vector of weights for each pigment
#' @param S_dvChl Optional vector of divinyl chlorophyll concentrations for Prochlorococcus
#' @param chemtax_style uncapped column-mean weighting with no chl_a 
#'          override (default = FALSE)
#'
#' @return A list containing the following elements:
#'   \item{F matrix}{The normalized F matrix of pigment ratios}
#'   \item{RMSE}{Root mean square error of the fit}
#'   \item{condition number}{Condition number of Fn %*% t(S)}
#'   \item{Class abundances}{Data frame of phytoplankton class abundances}
#'   \item{Figure}{Plot of the results}
#'   \item{MAE}{Mean absolute error for each pigment}
#'   \item{Error}{Residual error matrix}
#'
#' @examples
#'  Fmat <- as.matrix(phytoclass::Fm)
#'  S <- as.matrix(phytoclass::Sm)
#'  S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
#'  S_Chl <- S[, ncol(S)]
#'  # Run NNLS_MF_Final
#'  result <- phytoclass:::NNLS_MF_Final(Fmat, S, S_Chl, S_weights)
#' 
#' 
NNLS_MF_Final <- function(Fn, S, S_Chl, S_weights, S_dvChl = NULL,
                          chemtax_style = FALSE) {
  
  check_pro <- any(tolower(colnames(Fn)) %in% c("dvchl", "dvchla", "chlvp"))

  if (check_pro && chemtax_style) {
    warning(
      "chemtax_style not yet supported for prochloro samples; ",
      "falling back to phytoclass normalisation."
    )
    chemtax_style <- FALSE
  }

  if (check_pro) {
    F_norm <- Prochloro_Normalise_F(Fn)
    Fn <- F_norm[[1]] * F_norm[[2]]
  } else {
    F_norm <- Normalise_F(Fn, chemtax_style = chemtax_style)
    Fn <- if (chemtax_style) F_norm[[1]] else F_norm[[1]] * F_norm[[2]]
  }

  equality_sum <-
    if (chemtax_style) {
      1
    } else if (check_pro) {
      NULL
    } else {
      S[, ncol(S)]
    }

  C_new2 <- nnls_lsei(Fn, S, S_weights, equality_sum = equality_sum)

  if (check_pro) {
    n_cls                        <- ncol(C_new2)
    new_pro                      <- S_dvChl
    new_pro[!is.finite(new_pro)] <- 0
    new_pro                      <- pmin(pmax(new_pro, 0), S_Chl)
    C_nonpro                     <- C_new2[, -n_cls, drop = FALSE]
    Cn_nonpro_sum                <- rowSums(C_nonpro)
    Cn_nonpro_sum                <- ifelse(Cn_nonpro_sum == 0, 1, Cn_nonpro_sum)
    Cn_nonpro                    <- (C_nonpro / Cn_nonpro_sum) * (S_Chl - new_pro)
    Cn2                          <- cbind(Cn_nonpro, new_pro)
    colnames(Cn2)                <- rownames(Fn)
    
  } else if (chemtax_style) {
    F_chla              <- Fn[, ncol(Fn)]
    S_chla              <- S[, ncol(S)]
    S_chla[S_chla == 0] <- 1
    Cn2                 <- sweep(C_new2, 2, F_chla, "*")
    Cn2                 <- sweep(Cn2, 1, S_chla, "/")
    Cn2                 <- sweep(Cn2, 1, S_Chl, "*")
    Cn_sum              <- rowSums(Cn2)
    Cn_sum              <- ifelse(Cn_sum == 0, 1, Cn_sum)
    Cn2                 <- sweep(Cn2 / Cn_sum, 1, S_Chl, "*")
    colnames(Cn2)       <- rownames(Fn)
    
  } else {
    Cn.s2         <- rowSums(C_new2)
    Cn.s2         <- ifelse(Cn.s2 == 0, 1, Cn.s2)
    Cn2           <- (C_new2 / Cn.s2) * S_Chl
    colnames(Cn2) <- rownames(Fn)
  }

  colnames(Fn)  <- colnames(S)
  Cn2           <- as.data.frame(Cn2)
  rownames(Cn2) <- rownames(S)

  S_residual <- S - (C_new2 %*% Fn)
  S_rmse     <- sqrt(mean(S_residual^2))
  S_mae      <- colMeans(abs(S_residual))
  cd         <- kappa(Fn %*% t(S))
  plt        <- phyto_figure(Cn2)

  list(
    `F matrix`         = Fn,
    RMSE               = S_rmse,
    `condition number` = cd,
    `Class abundances` = Cn2,
    Figure             = plt,
    MAE                = S_mae,
    Error              = S_residual
  )
}
# ============================================================================ #
# ---- old versions ---- #
# ============================================================================ #
#NNLS_MF_Final <- function(Fn, S, S_Chl, S_weights, S_dvChl = NULL) {
#   check_pro <- any(tolower(colnames(Fn)) %in% c("dvchl", "dvchla", "chlvp"))
  
#   # normalize F matrix
#   if (check_pro) {
#     F_norm <- Prochloro_Normalise_F(Fn)
#   } else {
#     F_norm <- Normalise_F(Fn)
#   }
#   Fn     <- F_norm[[1]] * F_norm[[2]]
  
#   Fn_wt_err <- t(Weight_error(Fn, S_weights))
#   S_wt_err  <- t(Weight_error(S, S_weights))
  
#   b       <- crossprod(Fn_wt_err, S_wt_err) # right hand side of linear eq
#   Fn_prod <- crossprod(Fn_wt_err) # positive definite matrix with coefficients 
  
#   # ---- calc NNLS ---- #
#   C_new2  <- 
#     RcppML::nnls(
#       Fn_prod,
#       b,
#       cd_maxit = 1000, 
#       cd_tol   = 1e-10
#     )
#   C_new2 <- t(C_new2)
  
#   C_new2        <- as.matrix(C_new2)
#   Cn.s2         <- rowSums(C_new2)
#   Cn.s2         <- ifelse(Cn.s2 == 0, 1, Cn.s2)
#   Cn2           <- C_new2 / Cn.s2
#   Cn2           <- as.matrix(Cn2)
#   Cn2           <- Cn2 * S_Chl
#   colnames(Cn2) <- rownames(Fn)
#   colnames(Fn)  <- colnames(S)
#   Cn2           <- as.data.frame(Cn2)
#   rownames(Cn2) <- rownames(S)
  
#   if (check_pro) {
#     # determine final class abundance for prochloro
#     Fn_wt_err2 <- t(Weight_error(Fn[nrow(Fn), -ncol(Fn)], S_weights[-length(S_weights)]))
#     S_wt_err2  <- t(Weight_error(S[, -ncol(S)], S_weights[-length(S_weights)]))
    
#     Pb       <- crossprod(Fn_wt_err2, S_wt_err2)
#     Fn_prod2 <- crossprod(Fn_wt_err2)
    
#     PC_new2 <-
#       RcppML::nnls(
#         Fn_prod2,
#         Pb,
#         cd_maxit = 1000,
#         cd_tol   = 1e-10
#       )
#     PC_new2 <- t(PC_new2)
#     PCn.s2  <- rowSums(PC_new2)
#     PCn2    <- PC_new2 / PCn.s2
#     PCn2    <- PCn2 * S_dvChl

#     Cn2[, ncol(Cn2)] <- as.vector(PCn2)
#   }
  
#   # ---- calculate error terms ---- #
#   S_residual <- S - (C_new2 %*% Fn)       # residual error
#   S_rmse     <- sqrt(mean(S_residual^2))  # RMSE
#   S_mae      <- colMeans(abs(S_residual)) # MAE
  
#   # ---- condition number ---- #
#   cd <- kappa(Fn %*% t(S))
  
#   # ---- plot final results ---- #
#   plt <- phyto_figure(Cn2)
  
#   return(list(
#     "F matrix"         = Fn,
#     "RMSE"             = S_rmse,
#     "condition number" = cd,
#     "Class abundances" = Cn2,
#     "Figure"           = plt,
#     "MAE"              = S_mae,
#     "Error"            = S_residual
#   )) 
# }

#' #' Final step for MF with prochlorococcus
#' #' @keywords internal
#' #'
#' #' @param Fn 
#' #' @param S 
#' #' @param S_Chl 
#' #' @param S_weights 
#' #' @param S_dvChl 
#' #'
#' #' @return
#' #'
#' #' @examples
#' Prochloro_NNLS_MF_Final <- function(Fn, S, S_Chl, S_weights, S_dvChl) {
#'   # normalize F matrix
#'   F_norm <- Prochloro_Normalise_F(Fn)
#'   Fn     <- F_norm[[1]] * F_norm[[2]]
#'   
#'   
#'   # determine final class abundance for non-prochloro
#'   Fn_wt_err <- t(Weight_error(Fn, S_weights))
#'   S_wt_err  <- t(Weight_error(S, S_weights))
#'   
#'   b       <- crossprod(Fn_wt_err, S_wt_err) # right hand side of linear eq
#'   Fn_prod <- crossprod(Fn_wt_err) # positive definite matrix with coefficients 
#'   
#'   C_new2 <- 
#'     RcppML::nnls(
#'       Fn_prod,
#'       b, 
#'       cd_maxit = 1000, 
#'       cd_tol = 1e-10
#'     )
#'   C_new2 <- t(C_new2)
#'   C_new2 <- as.matrix(C_new2)
#'   Cn.s2  <- rowSums(C_new2)
#'   Cn2    <- C_new2 / Cn.s2
#'   Cn2    <- as.matrix(Cn2)
#'   Cn2    <- Cn2 * S_Chl
#'   
#'   colnames(Cn2) <- rownames(Fn)
#'   colnames(Fn)  <- colnames(S)
#'   
#'   # determine final class abundance for prochloro
#'   Fn_wt_err2 <- t(Weight_error(Fn[nrow(Fn), -ncol(Fn)], S_weights[-length(S_weights)]))
#'   S_wt_err2  <- t(Weight_error(S[, -ncol(S)], S_weights[-length(S_weights)]))
#'   
#'   Pb       <- crossprod(Fn_wt_err2, S_wt_err2)
#'   Fn_prod2 <- crossprod(Fn_wt_err2)
#'   
#'   PC_new2 <- 
#'     RcppML::nnls(
#'       Fn_prod2,
#'       Pb,
#'       cd_maxit = 1000,
#'       cd_tol   = 1e-10
#'     )
#'   PC_new2 <- t(PC_new2)
#'   PCn.s2  <- rowSums(PC_new2)
#'   PCn2    <- PC_new2 / PCn.s2
#'   PCn2    <- PCn2 * S_dvChl
#'   
#'   Cn2[, ncol(Cn2)] <- PCn2
#'   Cn2 <- as.data.frame(Cn2)
#'   
#'   # ---- calculate error terms ---- #
#'   S_residual <- S - (C_new2 %*% Fn)       # residual error
#'   S_rmse     <- sqrt(mean(S_residual^2))  # RMSE
#'   S_mae      <- colMeans(abs(S_residual)) # MAE
#'   
#'   # ---- condition number ---- #
#'   cd <- kappa(Fn %*% t(S))
#'   
#'   # ---- plot final results ---- #
#'   plt <- phyto_figure(Cn2)
#'   
#'   # row.names(Cn2) <- row.names(G)
#'   
#'   return(
#'     list(
#'       "F matrix"         = Fn, 
#'       "RMSE"             = S_rmse, 
#'       "condition number" = cd,
#'       "Class abundances" = Cn2, 
#'       "Figure"           = plt, 
#'       "MAE"              = S_mae, 
#'       "Error"            = S_residual
#'     )
#'   )
#' }
