test_that("Condition_test returns a numeric mean condition number for valid input", {
  # Example setup based on simulated_annealing usage
  Fmat <- as.matrix(phytoclass::Fm)
  S <- as.matrix(phytoclass::Sm)
  min_max <- phytoclass::min_max
  min.val <- min_max[[3]]
  max.val <- min_max[[4]]
  # Use only non-chla columns for Fmat and S, as in simulated_annealing
  Fmat_sub <- Fmat[, -ncol(Fmat)]
  S_sub <- S[, -ncol(S)]
  # Run Condition_test
  cond <- phytoclass:::Condition_test(S_sub, Fmat_sub, min.val, max.val)
  expect_type(cond, "double")
  expect_true(cond > 0)
})
