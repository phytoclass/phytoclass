test_that("Steepest_Descent returns expected output for valid input", {
  # Setup based on simulated_annealing usage
  Fmat <- as.matrix(phytoclass::Fm)
  S <- as.matrix(phytoclass::Sm)
  S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
  place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
  num.loops <- 2
  # Run Steepest_Descent
  result <- phytoclass:::Steepest_Descent(Fmat, place, S, S_weights, num.loops)
  # Should return a list with F matrix, RMSE, and C matrix
  expect_type(result, "list")
  expect_true(is.matrix(result[[1]]))
  expect_true(is.numeric(result[[2]]))
  expect_true(is.matrix(result[[3]]))
})
