#' Weighting
#'
#' @param Fn 
#'
#' @return A list
#' @export
#'
#' @examples
Fac_F <- function(Fn){
  b <- crossprod(t(Weight_error(Fn)),t(Weight_error(S)))
  C_new2 <-t(nnls(crossprod(t(Weight_error(Fn))),b,cd_maxit = 100000,cd_tol =1e-6 ))
  C_new2 <- as.matrix(C_new2)
  Cn.s2 <- rowSums(C_new2)
  Cn2 <- C_new2/Cn.s2 #Row sums to one
  Cn2 <- as.matrix(Cn2)
  colnames(Cn2) <- rownames(F)
  error <- norm(Weight_error(S)-(C_new2%*%Weight_error(Fn)),'F')
  return(list(Fn,error,Cn2))
}

