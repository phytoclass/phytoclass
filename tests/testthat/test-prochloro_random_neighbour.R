test_that("Prochloro_Random_Neighbour returns expected output for valid input", {
  # Setup based on simulated_annealing_Prochloro usage
  Fmat <- as.matrix(phytoclass::Fp)
  S <- as.matrix(phytoclass::Sp)
  S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
  place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
  
  # Get min_max from package data
  min_max_mat <- phytoclass:::Default_min_max(phytoclass::min_max, Fmat[, -ncol(Fmat)])
  
  # Get bounds using Prochloro_Wrangling as in simulated_annealing_Prochloro
  f_c <- Fmat
  W0 <- phytoclass:::Prochloro_Wrangling(f_c, min_max_mat[[1]], min_max_mat[[2]])
  minF <- W0[[1]]
  maxF <- W0[[2]]
  
  # Extract chlv and chlvp from the matrix
  chlv <- f_c[, ncol(f_c)]        # Chl a (Tchla)
  chlvp <- f_c[, ncol(f_c) - 1]   # dvChl a
  
  Temp <- 0.5
  N <- place
  
  # Run Prochloro_Random_Neighbour
  result <- phytoclass:::Prochloro_Random_Neighbour(
    f_c, Temp, chlv, chlvp, N, place, S, S_weights, minF, maxF
  )
  
  # Should return a list with F matrix, RMSE, and C matrix (from NNLS_MF)
  expect_type(result, "list")
  expect_length(result, 3)
  expect_true(is.matrix(result[[1]]))  # F matrix
  expect_true(is.numeric(result[[2]]))  # RMSE
  expect_true(is.matrix(result[[3]]))  # C matrix
})
