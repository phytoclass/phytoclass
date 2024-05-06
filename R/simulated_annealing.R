#' This is the main phytoclass algorithm. It performs simulated annealing algorithm for S and F matrices. See the examples (Fm, Sm) for how to set up matrices, and the vignette for more detailed instructions. Different pigments and phytoplankton groups may be used.
#'
#' @param S   Sample data matrix – a matrix of pigment samples
#' @param Fmat   Pigment to Chl a matrix
#' @param user_defined_min_max data frame with some format as min_max built-in data
#' @param do_matrix_checks     This should only be set to TRUE when using the default values. This will remove pigment columns that have column sums of 0. Set to FALSE if using customised names for pigments and phytoplankton groups
#' @param niter Number of iterations (default is 500)
#' @param step  Step ratio used (default is 0.009)
#' @param weight.upper.bound Upper limit of the weights applied (default value is 30). 
#' @param verbose Logical value. Output error and temperature at each iteration. Default value of TRUE
#'
#' @return A list containing 
#' \enumerate{
#'  \item Fmat matrix
#'  \item RMSE (Root Mean Square Error)
#'  \item condition number
#'  \item Class abundances
#'  \item Figure (plot of results)
#'  \item MAE (Mean Absolute Error)
#'  \item Error
#'  }
#' @export
#'
#' @examples
#' # Using the built-in matrices Sm and Fm
#' set.seed(5326)
#' sa.example <- simulated_annealing(Sm, Fm, niter = 5)
#' sa.example$Figure
simulated_annealing <- function(S,
                                Fmat = NULL, 
                                user_defined_min_max = NULL,
                                do_matrix_checks = TRUE,
                                niter = 500,
                                step = 0.009,
                                weight.upper.bound = 30, 
                                verbose = TRUE){
  if (is.null(Fmat)) {
    Fmat <- phytoclass::Fm
  }

  if (is.data.frame(S)){
  char_cols <- sapply(S, is.character)
  S <- S[, !char_cols]}

  if (do_matrix_checks) {
    L <- Matrix_checks(as.matrix(S), as.matrix(Fmat))
    S <- as.matrix(L[[1]])
    Fmat <- as.matrix(L[[2]])
  }
  
  S_Chl <- S[, ncol(S)]
  S <- Normalise_S(S)
  cm <- Bounded_weights(S, weight.upper.bound)
  place <- which(Fmat[,1:ncol(Fmat)-1] > 0)
  
  if (is.null(user_defined_min_max)) {
    K <- Default_min_max(phytoclass::min_max, Fmat[,1:ncol(Fmat)-1], place)
    min.val <- K[[1]]
    max.val <- K[[2]]
  }
  else {
    K <- Default_min_max(user_defined_min_max,Fmat[, 1:ncol(Fmat) - 1], place)
  min.val <- K[[1]]
  max.val <- K[[2]]
    # if (length(min.val) != length(place)) {
    #   message(paste0("\nNumber of rows for user_defined_min_max = ", 
    #                  length(min.val)))
    #   message(paste0("Length of place = ", length(place)))
    #   stop("\nThese do not match.")
    # }
  }
  
  condition.test <- Condition_test(S[,1:ncol(S)-1], Fmat[,1:ncol(Fmat)-1], min.val, max.val)
  if (verbose) {
   message(paste0(
     "\nCondition number = ", round(condition.test),
     "\n\n"
   ))
 }
  
  if (condition.test > 10^5) {
    stop("Condition number of S matrix greater than 100 000\n")
  }
  
  Fi <- ifelse(Fmat > 0, 1, 0)
  SE <- vectorise(Fi)
  nc <- NNLS_MF(Fi, S, cm)
  s_b <- s_c <- s_n <- nc[[1]]
  f_b <- f_c <- f_n <- nc[[2]]
  
  for (k in 1:niter) {
    Temp <- (1 - step)^(k)
    chlv <- Wrangling(s_c, min.val, max.val)[[4]]
    new_neighbour <- Random_neighbour2(s_c, Temp, chlv, s_c, 
                                       place, S, cm, min.val, max.val)
    D <- list()
    if (k > niter-20){
      for (i in 1:300){
        chlv <- Wrangling(s_c, min.val, max.val)[[4]]
        D[[length(D)+1]] <- Random_neighbour2(s_c, Temp,
                                              chlv, s_c, place, S, cm, min.val, max.val)
      }
    }
    else{
      for (i in 1:120){
        chlv <- Wrangling(s_c, min.val, max.val)[[4]]
        D[[length(D)+1]] <- Random_neighbour2(s_c, Temp,
                                              chlv, s_c, place, S, cm, min.val, max.val)
      }
    }
    Dn <- list()
    for (i in D) {
      Dn[[length(Dn) + 1]] <- i[[2]]
    }
    nk <- which.min(Dn)
    new_neighbour <- D[[nk]]
    if (Temp > 0.3) {
      new_neighbour <- SAALS(new_neighbour[[1]], min.val, 
                             max.val, place, S, cm, num.loops = 10)
    }
    else {
      new_neighbour <- SAALS(new_neighbour[[1]], min.val, 
                             max.val, place, S, cm, num.loops = 2)
    }
    minF <- Wrangling(new_neighbour[[1]], min.val, max.val)[[1]]
    maxF <- Wrangling(new_neighbour[[1]], min.val, max.val)[[2]]
    s_n <- new_neighbour[[1]]
    f_n <- new_neighbour[[2]]
    loop <- 1
    d <- which(vectorise(s_n[, 1:(ncol(s_n) - 1)]) < minF | 
                 vectorise(s_n[, 1:(ncol(s_n) - 1)]) > maxF)
    while (length(d) > 0) {
        if (k > niter-20){
          N <- place[d]
          for (i in 1:300){
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
      else{
        N <- place[d]
        for (i in 1:120){
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
      
    }
    
    A = target(f_n)/target(f_c)
    diff <- f_n - f_c
    if (f_n < f_c || exp(-(f_n - f_c)) < stats::runif(1, 
                                                      0, 1)) {
      s_c <- s_n
      f_c <- f_n
    }
    
    if (verbose) {
      message(paste("Current error: ", round(f_c, 4)))
      message(paste("Neighbour's error: ", round(f_n, 4)))
      message(paste("Temperature (%): ", round(Temp * 100, 2)))
      message(" ")
    }
    
    if (f_n < f_b) {
      s_b <- s_n
      f_b <- f_n
    }
  }
  
  res <- list(s_b, f_b)
  A <- res[[1]]
  
  final.results <- NNLS_MF_Final(A, S, S_Chl, cm)
  return(final.results)
}
