#' Perform simulated annealing algorithm for samples with divinyl chlorophyll and prochlorococcus. Divinyl chlorophyll must be the final column of both S and F matrices, with chlorophyll a the 2nd to last column. See how the example Sp and Fp matrices are organised.
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
#' # Using the built-in matrices Sp and Fp.
#' set.seed(5326)
#' sa.example <- simulated_annealing_Prochloro(Sp, Fp, niter = 1)
#' sa.example$Figure
simulated_annealing_Prochloro <- function (S, 
                                           Fmat = NULL,
                                           user_defined_min_max = NULL,
                                           do_matrix_checks = TRUE, 
                                           niter = 500, 
                                           step = 0.009, 
                                           weight.upper.bound = 30, 
                                           verbose = TRUE) 
{
  if (is.null(Fmat)) {
    Fmat <- phytoclass::Fp
  }
  if (is.data.frame(S)) {
    char_cols <- sapply(S, is.character)
    S <- S[, !char_cols]
  }
  if (do_matrix_checks) {
    L <- Matrix_checks(as.matrix(S), as.matrix(Fmat))
    S <- as.matrix(L[[1]])
    Fmat <- as.matrix(L[[2]])
  }
  S_Chl <- S[, ncol(S)]
  S_dvChl <- S[, ncol(S) - 1]
  S <- Normalise_S(S)
  cm <- Bounded_weights(S, weight.upper.bound)
  place <- which(Fmat[, 1:ncol(Fmat) - 1] > 0)
  if (is.null(user_defined_min_max)) {
    K <- Default_min_max(phytoclass::min_max, 
                                      Fmat[, 1:ncol(Fmat) - 1], place)
    min.val <- K[[1]]
    max.val <- K[[2]]
  }
  else {
    K <- Default_min_max(user_defined_min_max, 
                                      Fmat[, 1:ncol(Fmat) - 1], place)
    min.val <- K[[1]]
    max.val <- K[[2]]
  }
  condition.test <- Condition_test(S[, 1:ncol(S) - 
                                                    1], Fmat[, 1:ncol(Fmat) - 1], min.val, max.val)
  if (verbose) {
    message(paste0("\nCondition number = ", round(condition.test), 
                   "\n\n"))
  }
  if (condition.test > 10^5) {
    print("Abort process: condition number of S matrix greater than 100 000\n")
  }
  Fi <- ifelse(Fmat > 0, 1, 0)
  SE <- vectorise(Fi)
  nc <- NNLS_MF(Fi, S, cm)
  s_b <- s_c <- s_n <- nc[[1]]
  f_b <- f_c <- f_n <- nc[[2]]
  for (k in 1:niter) {
    Temp <- (1 - step)^(k)
    chlv <- Prochloro_Wrangling(s_c, min.val, max.val)[[4]]
    chlvp <- Prochloro_Wrangling(s_c, min.val, max.val)[[5]]
    new_neighbour <- Prochloro_Random_Neighbour_2(s_c, Temp, 
                                                  chlv, s_c, place, S, cm, min.val, max.val, chlvp)
    D <- list()
    if (k > niter - 20) {
      for (i in 1:300) {
        chlv <- Prochloro_Wrangling(s_c, min.val, max.val)[[4]]
        chlvp <- Prochloro_Wrangling(s_c, min.val, max.val)[[5]]
        D[[length(D) + 1]] <- Prochloro_Random_Neighbour_2(s_c, 
                                                           Temp, chlv, s_c, place, S, cm, min.val, max.val, 
                                                           chlvp)
      }
    }
    else {
      for (i in 1:120) {
        chlv <- Prochloro_Wrangling(s_c, min.val, max.val)[[4]]
        chlvp <- Prochloro_Wrangling(s_c, min.val, max.val)[[5]]
        D[[length(D) + 1]] <- Prochloro_Random_Neighbour_2(s_c, 
                                                           Temp, chlv, s_c, place, S, cm, min.val, max.val, 
                                                           chlvp)
      }
    }
    Dn <- list()
    for (i in D) {
      Dn[[length(Dn) + 1]] <- i[[2]]
    }
    nk <- which.min(Dn)
    new_neighbour <- D[[nk]]
    if (Temp > 0.3) {
      new_neighbour <- SAALS(new_neighbour[[1]], 
                                          min.val, max.val, place, S, cm, num.loops = 10)
    }
    else {
      new_neighbour <- SAALS(new_neighbour[[1]], 
                                          min.val, max.val, place, S, cm, num.loops = 2)
    }
    minF <- Prochloro_Wrangling(new_neighbour[[1]], min.val, 
                                max.val)[[1]]
    maxF <- Prochloro_Wrangling(new_neighbour[[1]], min.val, 
                                max.val)[[2]]
    s_n <- new_neighbour[[1]]
    f_n <- new_neighbour[[2]]
    loop <- 1
    d <- which(vectorise(s_n[, 1:(ncol(s_n) - 
                                                 2)]) < minF[1:length(minF) - 1] | vectorise(s_n[, 
                                                                                                              1:(ncol(s_n) - 2)]) > maxF[1:length(maxF) - 1])
    while (length(d) > 0) {
      if (k > niter - 20) {
        N <- place[d]
        for (i in 1:300) {
          chlv <- Prochloro_Wrangling(s_n, min.val, max.val)[[4]]
          chlvp <- Prochloro_Wrangling(s_c, min.val, 
                                       max.val)[[5]]
          D[[length(D) + 1]] <- Prochloro_Random_Neighbour(s_n, 
                                                           Temp, chlv, s_n, N, place, S, cm, min.val, 
                                                           max.val)
        }
        Dn <- list()
        for (i in D) {
          Dn[[length(Dn) + 1]] <- i[[2]]
        }
        nk <- which.min(Dn)
        new_neighbour <- D[[nk]]
        s_n <- new_neighbour[[1]]
        f_n <- new_neighbour[[2]]
        minF <- Prochloro_Wrangling(new_neighbour[[1]], 
                                    min.val, max.val)[[1]]
        maxF <- Prochloro_Wrangling(new_neighbour[[1]], 
                                    min.val, max.val)[[2]]
        d <- which(vectorise(s_n[, 1:(ncol(s_n) - 
                                                     2)]) < minF[1:length(minF) - 1] | vectorise(s_n[, 
                                                                                                                  1:(ncol(s_n) - 2)]) > maxF[1:length(maxF) - 
                                                                                                                                               1])
      }
      else {
        N <- place[d]
        for (i in 1:120) {
          chlv <- Prochloro_Wrangling(s_n, min.val, max.val)[[4]]
          chlvp <- Prochloro_Wrangling(s_c, min.val, 
                                       max.val)[[5]]
          D[[length(D) + 1]] <- Prochloro_Random_Neighbour(s_n, 
                                                           Temp, chlv, s_n, N, place, S, cm, min.val, 
                                                           max.val)
        }
        Dn <- list()
        for (i in D) {
          Dn[[length(Dn) + 1]] <- i[[2]]
        }
        nk <- which.min(Dn)
        new_neighbour <- D[[nk]]
        s_n <- new_neighbour[[1]]
        f_n <- new_neighbour[[2]]
        minF <- Prochloro_Wrangling(new_neighbour[[1]], 
                                    min.val, max.val)[[1]]
        maxF <- Prochloro_Wrangling(new_neighbour[[1]], 
                                    min.val, max.val)[[2]]
        d <- which(vectorise(s_n[, 1:(ncol(s_n) - 
                                                     2)]) < minF[1:length(minF) - 1] | vectorise(s_n[, 
                                                                                                                  1:(ncol(s_n) - 2)]) > maxF[1:length(maxF) - 
                                                                                                                                               1])
      }
    }
    diff <- f_n - f_c
    if (f_n < f_c || exp(-(f_n - f_c)) < stats::runif(1, 
                                                      0, 1)) {
      s_c <- s_n
      f_c <- f_n
    }
    if (verbose) {
      message(paste("Current error: ", round(f_c, 4)))
      message(paste("Neighbour's error: ", round(f_n, 4)))
      message(paste("Temperature (%): ", round(Temp * 100, 
                                               2)))
      message(" ")
    }
    if (f_n < f_b) {
      s_b <- s_n
      f_b <- f_n
    }
  }
  res <- list(s_b, f_b)
  A <- res[[1]]
  final.results <- Prochloro_NNLS_MF_Final(A, S, S_Chl, cm, 
                                           S_dvChl)
  return(final.results)
}