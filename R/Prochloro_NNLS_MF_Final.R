#' Final step for MF with prochlorococcus
#' @keywords internal
#'
#' @param Fn 
#' @param S 
#' @param S_Chl 
#' @param cm 
#'
#' @return
#'
#' @examples
Prochloro_NNLS_MF_Final <- function (Fn, S, S_Chl, cm, S_dvChl) 
{
  F.sum <- Prochloro_Normalise_F(Fn)[[2]]
  Fn <- Prochloro_Normalise_F(Fn)[[1]]
  Fn <- Fn * F.sum
  b <- crossprod(t(Weight_error(Fn, cm)), t(Weight_error(S, 
                                                         cm)))
  C_new2 <- t(RcppML::nnls(crossprod(t(Weight_error(Fn, cm))), 
                           b, cd_maxit = 1000, cd_tol = 1e-10))
  C_new2 <- as.matrix(C_new2)
  Cn.s2 <- rowSums(C_new2)
  Cn2 <- C_new2/Cn.s2
  Cn2 <- as.matrix(Cn2)
  Cn2 <- Cn2 * S_Chl
  colnames(Cn2) <- rownames(Fn)
  colnames(Fn) <- colnames(S)
  Pb <- crossprod(t(Weight_error(Fn[nrow(Fn),1:ncol(Fn)-1], cm[1:length(cm)-1])), t(Weight_error(S[,1:ncol(S)-1],                                                                                                                       cm[1:length(cm)-1])))
  PC_new2 <- t(RcppML::nnls(crossprod(t(Weight_error(Fn[nrow(Fn),1:ncol(Fn)-1], cm[1:length(cm)-1]))), 
                            Pb, cd_maxit = 1000, cd_tol = 1e-10))
  PCn.s2 <- rowSums(PC_new2)
  PCn2 <- PC_new2/PCn.s2
  PCn2 <- PCn2 * S_dvChl
  Cn2[,ncol(Cn2)] <- PCn2
  error <- Metrics::rmse(S, (C_new2 %*% Fn))
  k <- Cn2
  k <- as.data.frame(k)
  k[, ncol(k) + 1] <- 1:nrow(k)
  cn <- colnames(k)
  cn <- cn[1:ncol(k) - 1]
  cd <- kappa(Fn %*% t(S))
  vals <- NULL
  .data <- NULL
  PLE <- tidyr::pivot_longer(data = k, cols = cn, names_to = "names", 
                             values_to = "vals")
  colorBlindGrey8 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#009E73", 
                       "#001E73", "#013E73")
  gr <- colnames(PLE)[[1]]
  n <- ggplot2::ggplot(PLE, ggplot2::aes(x = .data[[gr]], y = vals, 
                                         fill = names)) + ggplot2::geom_area() + ggplot2::scale_color_manual(values = colorBlindGrey8) + 
    ggplot2::scale_fill_manual(values = colorBlindGrey8) + 
    ggplot2::xlab("Sample number") + ggplot2::ylab("Chl a concentrations") + 
    ggplot2::theme_bw()
  G <- S - (C_new2 %*% Fn)
  Cn2 <- as.data.frame(Cn2)
  row.names(Cn2) <- row.names(G)
  gs <- colMeans(abs((C_new2%*%Fn) - S))
  return(list(`F matrix` = Fn, RMSE = error, `condition number` = cd, 
              `Class abundances` = Cn2, Figure = n, MAE = gs, Error = G))
}



