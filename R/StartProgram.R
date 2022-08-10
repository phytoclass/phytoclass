StartProgram <- function(S){
  Fmax <- ifelse(F>0,1,0)
  A <- simulated_annealing2(Fmax,niter = 450,0.009)[[1]]
  E <- Fac_F_Final(A)[[1]]
  return(E)
}