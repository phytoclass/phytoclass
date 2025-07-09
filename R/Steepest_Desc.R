#' Stand-alone version of steepest descent algorithm. This is similar to the CHEMTAX steepest descent algorithm. It is not required to use this function, and as results are not bound by minimum and maximum, results may be unrealistic.
#'
#' @param S   Sample data matrix – a matrix of pigment samples
#' @param Fmat   Pigment to Chl a matrix
#' @param num.loops Number of loops/iterations to perform (no default)
#'
#' @return A list containing 
#' \enumerate{
#'  \item The F matrix (pigment: Chl *a*) ratios
#'  \item RMSE (Root Mean Square Error)
#'  \item Condition number
#'  \item class abundances
#'  \item Figure (plot of results)
#'  \item MAE (Mean Absolute Error)
#'  }
#' @export
#'
#' @examples
#' MC <- Matrix_checks(Sm,Fm)
#' Snew <- MC$Snew
#' Fnew <- MC$Fnew
#' SDRes <- Steepest_Desc(Fnew,Snew, num.loops = 20)
#' 
Steepest_Desc <- function(Fmat, S, num.loops) {
  S_Chl <- S[, ncol(S)]
  cm    <- Bounded_weights(S, 30)
  place <- which(Fmat[, -ncol(Fmat)] > 0)

  # loop      <- 1
  F.new     <- NNLS_MF(Fmat, S, cm)
  F.initial <- F.new
  
  for (i in 1:num.loops) {
    F.new  <- Minimise_elements_comb(F.initial[[1]], place, S, cm, c1_num = 3)
    # loop   <- loop + 1
    loop <- 1
    while (F.new[[2]] > F.initial[[2]]) {
      
      if (loop <= 5) {
        F.new <- Minimise_elements_comb(F.initial[[1]], place, S, cm, c1_num = 3)
      } else if (loop < 10) {
        F.new <- Minimise_elements_comb(F.initial[[1]], place, S, cm, c1_num = 1)
      } else if (loop <= 100) {
        # If it doesn't work the first time, it randomises at a lower rate
        F.new <- Minimise_elements_comb(F.initial[[1]], place, S, cm, c1_num = 2)
      } else if (loop > 100) {
        # it will continue for 100 iterations, and then stop
        break
      }
      loop <- loop + 1
    }
    F.initial <- F.new
  }
  F.new <- NNLS_MF_Final(as.matrix(F.new[[1]]), S, S_Chl, cm)
  return(F.new)
}



