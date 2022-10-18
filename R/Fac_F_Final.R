#' This is the final fit, and multiplies the final value by Chl concentration so 
#' they're in units of Chl a biomass
#'
#' @param Fn 
#' @param S
#' @param S_Chl
#'
#' @return
#' @export
#'
#' @examples
Fac_F_Final <- function(Fn, S, S_Chl){
  F.sum <- Wrangle_F(Fn)[[2]]
  Fn <- Wrangle_F(Fn)[[1]]
  Fn <- Fn * F.sum
  b <- crossprod(t(Fn),t(S))
  C_new2 <-t(RcppML::nnls(crossprod(t(Fn)),b,cd_maxit = 1000000000,cd_tol =1e-18 ))
  C_new2 <- as.matrix(C_new2)
  Cn.s2 <- rowSums(C_new2)
  Cn2 <- C_new2/Cn.s2
  Cn2 <- as.matrix(Cn2)
  Cn2 <- Cn2 * S_Chl
  colnames(Cn2) <- rownames(Fn)
  error <- norm(S-(C_new2%*%Fn),'F')
  return(list(Fn,error,Cn2))
}
