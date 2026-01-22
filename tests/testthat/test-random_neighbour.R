test_that("Random_neighbour returns expected output for valid input", {
  # Setup based on simulated_annealing usage
  Fmat <- as.matrix(phytoclass::Fm)
  S <- as.matrix(phytoclass::Sm)
  S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
  place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
  min_max <- phytoclass::min_max
  minF <- min_max[[3]][seq_along(place)]
  maxF <- min_max[[4]][seq_along(place)]
  chlv <- rep(1, nrow(Fmat)) # typical usage in simulated_annealing
  Temp <- 0.5
  N <- place
  # Run Random_neighbour
  result <- phytoclass:::Random_neighbour(Fmat, Temp, chlv, N, place, S, S_weights, minF, maxF)
  # Should return a list with F matrix, RMSE, and C matrix
  expect_type(result, "list")
  expect_true(is.matrix(result[[1]]))
  expect_true(is.numeric(result[[2]]))
  expect_true(is.matrix(result[[3]]))
})
