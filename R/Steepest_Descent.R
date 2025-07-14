#' Performs the steepest descent algorithm for a set number of iterations
#' 
#' @keywords internal
#' 
#' @param Fmat xx
#' @param place xx
#' @param S   xx
#' @param S_weights   xx
#' @param num.loops   xx
#'
#' @return
#'
#' @examples
Steepest_Descent <- function(Fmat, place, S, S_weights, num.loops) {
  F_new     <- NNLS_MF(Fmat, S, S_weights)
  F_initial <- F_new
  
  for (i in 1:num.loops) { # should always be small. It would be nice to allow the
    F_new <- Minimise_elements_comb(F_initial[[1]], place, S, S_weights, c1_num = 3)

    loop <- 1
    while (F_new[[2]] > F_initial[[2]]) {

      if (loop <= 5) {
        F_new <- Minimise_elements_comb(F_initial[[1]], place, S, S_weights, c1_num = 3)
      } else if (loop < 10) {
        F_new <- Minimise_elements_comb(F_initial[[1]], place, S, S_weights, c1_num = 1)
      } else if (loop <= 100) {
        # If it doesn't work the first time, it randomises at a lower rate
        F_new <- Minimise_elements_comb(F_initial[[1]], place, S, S_weights, c1_num = 2)
      } else if (loop > 100) {
        # it will continue for 100 iterations, and then stop
        break
      }
      
      loop <- loop + 1
    }
    F_initial <- F_new
  }
  return(F_new)
}
