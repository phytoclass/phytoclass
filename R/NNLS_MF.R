#' Performs the non-negative matrix factorisation for given phytoplankton 
#' pigments and pigment ratios, to attain an estimate of phytoplankton 
#' class abundances.
#'
#' @param Fn   xx
#' @param S    xx
#' @param cm  xx
#'
#' @return A list
#' @export
#'
#' @examples
#'
NNLS_MF <- function(Fn, S, cm=NULL){
    if (is.null(cm)) {
    cm <- as.vector(rep(1,ncol(S)))
  }


  b <- crossprod(t(Weight_error(Fn, cm)),t(Weight_error(S, cm)))
  C_new2 <-t(RcppML::nnls(crossprod(t(Weight_error(Fn, cm))),
                          b, 
                          cd_maxit = 1000,cd_tol =1e-8 ))
  C_new2 <- as.matrix(C_new2)
  Cn.s2 <- rowSums(C_new2)
  Cn2 <- C_new2/Cn.s2 #Row sums to one
  Cn2 <- as.matrix(Cn2)
  colnames(Cn2) <- rownames(Fn)
  error <- Metrics::rmse((S),C_new2%*%(Fn))
  return(list("F matrix " = Fn, "RMSE" = error,"C matrix" = Cn2))
}

