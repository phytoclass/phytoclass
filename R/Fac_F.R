#' Weighting
#'
#' @param Fn 
#' @param S 
#' @param cm
#'
#' @return A list
#' @export
#'
#' @examples
#'
Fac_F <- function(Fn, S, cm){
  b <- crossprod(t(Weight_error(Fn, cm)),t(Weight_error(S, cm)))
  C_new2 <-t(RcppML::nnls(crossprod(t(Weight_error(Fn, cm))), b, 
                          cd_maxit = 100000,cd_tol =1e-6 ))
  C_new2 <- as.matrix(C_new2)
  Cn.s2 <- rowSums(C_new2)
  Cn2 <- C_new2/Cn.s2 #Row sums to one
  Cn2 <- as.matrix(Cn2)
  colnames(Cn2) <- rownames(F)
  error <- norm(Weight_error(S, cm) - (C_new2%*%Weight_error(Fn, )), 'F')
  return(list(Fn, error, Cn2))
}

