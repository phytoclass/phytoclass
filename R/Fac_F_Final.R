#' This is the final fit, and multiplies the final value by Chl concentration so 
#' they're in units of Chl a biomass
#'
#' @param Fn 
#' @param S
#' @param S_Chl
#' @param cm
#'
#' @return
#' @export
#'
#' @examples
Fac_F_Final <- function(Fn, S, S_Chl, cm){
  F.sum <- Wrangle_F(Fn)[[2]]
  Fn <- Wrangle_F(Fn)[[1]]
  Fn <- Fn * F.sum
  b <- crossprod(t(Weight_error(Fn, cm)),t(Weight_error(S, cm)))
  C_new2 <-t(RcppML::nnls(crossprod(t(Weight_error(Fn, cm))), b,cd_maxit = 1000,cd_tol =1e-8 ))
  C_new2 <- as.matrix(C_new2)
  Cn.s2 <- rowSums(C_new2)
  Cn2 <- C_new2/Cn.s2
  Cn2 <- as.matrix(Cn2)
  Cn2 <- Cn2 * S_Chl
  colnames(Cn2) <- rownames(Fn)
  error <- Metrics::rmse(S,(C_new2%*%Fn))
  k <- Cn2
  k <- as.data.frame(k)
  k[,ncol(k)+1] <- 1:nrow(k)
  cn <- colnames(k)
  
  cn <- cn[1:ncol(k)-1]
  
  PLE <- tidyr::pivot_longer(data = k,cols = cn,names_to = 'names',values_to = 'vals' )
  colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#009E73","#001E73","#013E73")
  gr <- colnames(PLE)[[1]]
  n <- ggplot2::ggplot(PLE, aes(x=UQ(as.name(gr)), y=vals, fill=names)) +
    ggplot2::geom_area() +
    ggplot2::scale_color_manual(values=colorBlindGrey8)+
    ggplot2::scale_fill_manual(values=colorBlindGrey8)
  
  G <- S - (C_new2%*%Fn)
  gs <- colMeans(abs(G))/colSums(S)
  
  return(list(Fn,error,Cn2,n,gs,G))
}
