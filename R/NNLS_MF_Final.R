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
#' @param cm  xx
#'
#' @return
#'
#' @examples  
NNLS_MF_Final <- function(Fn, S, S_Chl, cm){
  F.sum <- Normalise_F(Fn)[[2]]
  Fn <- Normalise_F(Fn)[[1]]
  Fn <- Fn * F.sum
  b <- crossprod(t(Weight_error(Fn, cm)),t(Weight_error(S, cm)))
  
  C_new2 <-t(RcppML::nnls(crossprod(t(Weight_error(Fn, cm))), b,cd_maxit = 1000,cd_tol =1e-10 ))
  C_new2 <- as.matrix(C_new2)
  Cn.s2 <- rowSums(C_new2)
  Cn2 <- C_new2/Cn.s2
  Cn2 <- as.matrix(Cn2)
  Cn2 <- Cn2 * S_Chl
  colnames(Cn2) <- rownames(Fn)
  colnames(Fn) <- colnames(S)
  
  error <- Metrics::rmse(S,(C_new2%*%Fn))

  k <- Cn2
  k <- as.data.frame(k)
  k[,ncol(k)+1] <- 1:nrow(k)
  
  cn <- colnames(k)
  cn <- cn[1:ncol(k)-1]
  
  # condition number
  cd <- kappa(Fn %*% t(S))
  
  # NULL assignment to stop NOTE during the package "Check"
  #  -  no visible binding for global variable
  vals <-  NULL
  .data <- NULL
  PLE <- tidyr::pivot_longer(data = k, cols = cn, names_to = 'names', values_to = 'vals' )
  colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#009E73","#001E73","#013E73")
  gr <- colnames(PLE)[[1]]
  
  n <- ggplot2::ggplot(PLE, ggplot2::aes(x = .data[[gr]], y=vals, fill=names)) +
    ggplot2::geom_area() +
    ggplot2::scale_color_manual(values=colorBlindGrey8)+
    ggplot2::scale_fill_manual(values=colorBlindGrey8) +
    ggplot2::xlab("Sample number") +
    ggplot2::ylab("Chl a concentrations") +
    ggplot2::theme_bw()
  
  G <- S - (C_new2%*%Fn)
  gs <- colMeans(abs((C_new2%*%Fn) - S))
  
  return(list("F matrix" = Fn, 
              "RMSE"  = error,
              "condition number" = cd,
              "Class abundances" = Cn2,
              "Figure" = n, 
              "MAE" = gs, 
              "Error" = G))  
}
