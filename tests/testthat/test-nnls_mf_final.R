test_that("NNLS_MF_Final returns expected output for valid input", {
  # Setup based on simulated_annealing usage
  Fmat <- as.matrix(phytoclass::Fm)
  S <- as.matrix(phytoclass::Sm)
  S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
  S_Chl <- S[, ncol(S)]
  # Run NNLS_MF_Final
  result <- phytoclass:::NNLS_MF_Final(Fmat, S, S_Chl, S_weights)
  # Should return a list with expected named elements
  expect_type(result, "list")
  expect_true("F matrix" %in% names(result))
  expect_true("RMSE" %in% names(result))
  expect_true("Class abundances" %in% names(result))
  expect_true(is.matrix(result[["F matrix"]]))
  expect_true(is.numeric(result[["RMSE"]]))
  expect_true(is.data.frame(result[["Class abundances"]]))
})
