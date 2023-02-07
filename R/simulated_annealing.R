#' input = F matrix, number of iterations, and the step to use. 
#'
#' @param S
#' @param F
#' @param min_max
#' @param niter 
#' @param step  
#'
#' @return
#' @export
#'
#' @examples
simulated_annealing <- function(S, F, min_max, niter, step){
  
  L <- Homogenise_matrices(S, F)
  
  S <- as.matrix(L[[1]])
  S_Chl <- S[, ncol(S)]
  cm <- Bounded_weights(S)
  S <- Wrangle_S(S)
  
  F <- as.matrix(L[[2]])
  place <- which(F > 0)
  
  K <- minvalmaxval(min_max, F, place)
  min.val <- K[[1]]
  max.val <- K[[2]]
  Fi <- ifelse(F>0, 1, 0)
  
  
  SE <- vectorise(Fi)
  nc <- NNLS_MF(Fi, S, cm)
  
  s_b <- s_c <- s_n <- nc[[1]]  # sets initial values 
  f_b <- f_c <- f_n <- nc[[2]]

  for (k in 1:niter) #Set up loop (niter = number of iterations
  {   
    Temp <- (1 - step)^(k) ### Set temp to decline with each iteration
    #consider random neighbour
    chlv <- Wrangling(s_c, min.val, max.val)[[4]]
    new_neighbour <- Random_neighbour2(s_c, Temp, chlv, s_c, place, S, cm, min.val, max.val)
    
    # considers 50 random neighbours and picks the one with the lowest error
    D <- list()
    for (i in 1:70){
      chlv <- Wrangling(s_c, min.val, max.val)[[4]]
      D[[length(D)+1]] <- Random_neighbour2(s_c, Temp, chlv, s_c, place, S, cm, min.val, max.val)
    }
    Dn <- list()
    for (i in D){
      Dn[[length(Dn)+1]] <-  i[[2]]
    }
    
    nk <- which.min(Dn)
    
    new_neighbour <- D[[nk]]
    
    # Applies the steepst descent algorithm to the new neighbour
    
    if (Temp > .3) {
      new_neighbour <- SAALS(new_neighbour[[1]], min.val, max.val,
                             place, S, cm, num.loops = 10)
    }
    else {
      new_neighbour <- SAALS(new_neighbour[[1]], min.val, max.val, 
                              place, S, cm, num.loops = 2)
    }
    
    
    minF <- Wrangling(new_neighbour[[1]], min.val, max.val)[[1]] # mins and maxs
    maxF <- Wrangling(new_neighbour[[1]], min.val, max.val)[[2]]
    
    s_n <- new_neighbour[[1]]
    f_n <- new_neighbour[[2]]
    
    #Shows which are below / above the min and max values
    loop <- 1
    d <- which(vectorise(s_n[,1:(ncol(s_n)-1)])<minF | vectorise(s_n[,1:(ncol(s_n)-1)]) > maxF) 
    while( length(d) > 0){
      N <- place[d]
      for (i in 1:70){
        chlv <- Wrangling(s_n, min.val, max.val)[[4]]
        D[[length(D)+1]] <- Random_neighbour(s_n, Temp, chlv, s_n, N, place, S, cm, min.val, max.val)
      }
      Dn <- list()
      for (i in D){
        Dn[[length(Dn)+1]] <-  i[[2]]
      }
      
      nk <- which.min(Dn)
      
      new_neighbour <- D[[nk]]
      s_n <- new_neighbour[[1]]
      f_n <- new_neighbour[[2]]
      
      minF <- Wrangling(new_neighbour[[1]], min.val, max.val)[[1]]
      maxF <- Wrangling(new_neighbour[[1]], min.val, max.val)[[2]]
      
      d <- which(vectorise(s_n[,1:(ncol(s_n)-1)])<minF | vectorise(s_n[,1:(ncol(s_n)-1)]) > maxF) 
      
      
    }
    
    # Difference in error
    A = target(f_n)/target(f_c) 
    diff <- f_n - f_c
    if (f_n < f_c || exp(-(f_n - f_c)) < stats::runif(1, 0, 1)) {   # metropolis criterion we can switch off and on.
      s_c <- s_n
      f_c <- f_n
    }
    print(paste("Current: ", round(f_c,4))) # print if verbose
    print(paste("Neighbour: ", round(f_n,4))) # print if verbose
    print(paste("Temp: ",round(Temp,4)))
    print(" ")
    
    # update best state
    if (f_n < f_b) {
      s_b <- s_n
      f_b <- f_n       
    }
    
  }
  res <- list(s_b, f_b)
  A <- res[[1]]
 
  NNLS_MF_Final(A, S, S_Chl, cm)
  
}


