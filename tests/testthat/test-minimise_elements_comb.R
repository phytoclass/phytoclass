test_that("Minimise_elements_comb returns expected output for valid input", {
  # Setup based on Steepest_Descent usage
  Fmat <- as.matrix(phytoclass::Fm)
  S <- as.matrix(phytoclass::Sm)
  S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
  place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
  
  # Get F_initial from NNLS_MF as done in Steepest_Descent
  F_initial <- phytoclass::NNLS_MF(Fmat, S, S_weights)
  
  # Run Minimise_elements_comb with c1_num = 3 (as in Steepest_Descent)
  result <- phytoclass:::Minimise_elements_comb(F_initial[[1]], place, S, S_weights, c1_num = 3)
  
  # Should return a list with F matrix, RMSE, and C matrix
  expect_type(result, "list")
  expect_true(is.matrix(result[[1]]))
  expect_true(is.numeric(result[[2]]))
  expect_true(is.matrix(result[[3]]))
})

test_that("Minimise_elements_comb works with different c1_num values", {
  # Setup
  Fmat <- as.matrix(phytoclass::Fm)
  S <- as.matrix(phytoclass::Sm)
  S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
  place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
  F_initial <- phytoclass::NNLS_MF(Fmat, S, S_weights)
  
  # Test with c1_num = 1
  result_1 <- phytoclass:::Minimise_elements_comb(F_initial[[1]], place, S, S_weights, c1_num = 1)
  expect_type(result_1, "list")
  
  # Test with c1_num = 2
  result_2 <- phytoclass:::Minimise_elements_comb(F_initial[[1]], place, S, S_weights, c1_num = 2)
  expect_type(result_2, "list")
  
  # Test with c1_num = 3
  result_3 <- phytoclass:::Minimise_elements_comb(F_initial[[1]], place, S, S_weights, c1_num = 3)
  expect_type(result_3, "list")
})
