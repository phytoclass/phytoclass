test_that("Fac_F_RR returns expected output for valid input", {
  # Setup based on Minimise_elements_comb usage
  Fmat <- as.matrix(phytoclass::Fm)
  S <- as.matrix(phytoclass::Sm)
  cm <- as.numeric(phytoclass:::Bounded_weights(S))
  place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
  
  # Get F.new from Conduit as done in Minimise_elements_comb
  f <- phytoclass:::Conduit(Fmat, place, S, cm, c_num = 3)
  F.new <- f[[1]]
  
  # Set place1 as done in Minimise_elements_comb
  place1 <- place  # place1 = place when c1_num != 1
  
  # Run Fac_F_RR
  result <- phytoclass:::Fac_F_RR(F.new, vary = place, place = place1, S, cm, fac_rr = 3)
  
  # Should return a list with two elements
  expect_type(result, "list")
  expect_length(result, 2)
  # First element is updated F matrix (as a list from NNLS_MF)
  expect_type(result[[1]], "list")
  # Second element is vector of indices where improvements were found
  # Can be integer or numeric vector
  expect_true(is.vector(result[[2]]))
})

test_that("Fac_F_RR works with different fac_rr values", {
  # Setup
  Fmat <- as.matrix(phytoclass::Fm)
  S <- as.matrix(phytoclass::Sm)
  cm <- as.numeric(phytoclass:::Bounded_weights(S))
  place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
  
  # Get F.new from NNLS_MF
  F.new <- phytoclass::NNLS_MF(Fmat, S, cm)
  
  # Test with fac_rr = 1 (place1 = NULL when c1_num = 1)
  result_1 <- phytoclass:::Fac_F_RR(F.new, vary = place, place = NULL, S, cm, fac_rr = 1)
  expect_type(result_1, "list")
  expect_length(result_1, 2)
  
  # Test with fac_rr = 2 (place1 = place when c1_num != 1)
  result_2 <- phytoclass:::Fac_F_RR(F.new, vary = place, place = place, S, cm, fac_rr = 2)
  expect_type(result_2, "list")
  expect_length(result_2, 2)
  
  # Test with fac_rr = 3 (place1 = place when c1_num != 1)
  result_3 <- phytoclass:::Fac_F_RR(F.new, vary = place, place = place, S, cm, fac_rr = 3)
  expect_type(result_3, "list")
  expect_length(result_3, 2)
})
