#' Run the algorithm
#'
#' @param S 
#'
#' @return
#' @export
#'
#' @examples
StartProgram <- function(S, min.val, max.val){
  Fmax <- ifelse(F>0,1,0)
  A <- simulated_annealing2(Fmax, niter = 450,0.009, min.val, max.val)[[1]]
  E <- Fac_F_Final(A)[[1]]
  return(E)
}