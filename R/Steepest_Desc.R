#' Stand-alone version of steepest descent algorithm
#'
#' @param F xx
#' @param S xx
#' @param num.loops xx
#'
#' @return
#' @export
#'
#' @examples
Steepest_Desc <-  function (F, S, num.loops) 
{
  L <- Matrix_checks(S, F)
  S <- as.matrix(L[[1]])
  S_Chl <- S[, ncol(S)]
  cm <- Bounded_weights(S)
  S <- Normalise_S(S)
  F <- as.matrix(L[[2]])
  place <- which(F > 0)
  
  loop <- 1
  F.new <- NNLS_MF(F, S, cm)
  F.initial <- F.new
  for (i in 1:num.loops) {
    F.new <- Minimise_elements(F.initial[[1]], place, S, 
                               cm)
    loop = loop + 1
    loop_2 <- 1
    while (F.new[[2]] > F.initial[[2]]) {
      loop_2 = loop_2 + 1
      F.new <- Minimise_elements(F.initial[[1]], place, 
                                 S, cm)
      if (loop_2 > 5) {
        F.new <- Minimise_elements1(F.initial[[1]], place, 
                                    S, cm)
      }
      if (loop_2 > 10) {
        F.new <- Minimise_elements2(F.initial[[1]], place, 
                                    S, cm)
      }
      if (loop_2 > 100) {
        break
      }
    }
    F.initial <- F.new
  }
  F.new <- NNLS_MF_Final(F.new, S, S_Chl, cm)
  return(F.new)
}