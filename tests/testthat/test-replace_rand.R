test_that("Replace_Rand returns expected output for valid input", {
  # Setup based on Fac_F_RR usage
  Fmat <- as.matrix(phytoclass::Fm)
  S <- as.matrix(phytoclass::Sm)
  cm <- as.numeric(phytoclass:::Bounded_weights(S))
  
  # Get Fmat as a list from NNLS_MF (as used in Fac_F_RR)
  Fmat_list <- phytoclass::NNLS_MF(Fmat, S, cm)
  
  # Test with a single index
  i <- 1  # first non-zero element to modify
  min.scaler <- 0.99
  max.scaler <- 1.01
  
  # Run Replace_Rand
  result <- phytoclass:::Replace_Rand(Fmat_list, i, S, cm, min.scaler, max.scaler)
  
  # Should return a list with 4 elements (c() flattens the NNLS_MF list + logical)
  # The result is: list(F matrix, RMSE, C matrix, improved_logical)
  expect_type(result, "list")
  expect_length(result, 4)
  # First element is the F matrix
  expect_true(is.matrix(result[[1]]))
  # Second element is RMSE
  expect_true(is.numeric(result[[2]]))
  # Third element is C matrix
  expect_true(is.matrix(result[[3]]))
  # Fourth element is logical indicating if RMSE improved
  expect_type(result[[4]], "logical")
})

test_that("Replace_Rand works with different scaling factors", {
  # Setup
  Fmat <- as.matrix(phytoclass::Fm)
  S <- as.matrix(phytoclass::Sm)
  cm <- as.numeric(phytoclass:::Bounded_weights(S))
  Fmat_list <- phytoclass::NNLS_MF(Fmat, S, cm)
  i <- 1
  
  # Test with fac_rr = 1 scalers (0.99, 1.01)
  result_1 <- phytoclass:::Replace_Rand(Fmat_list, i, S, cm, 
                                        min.scaler = 0.99, max.scaler = 1.01)
  expect_type(result_1, "list")
  expect_length(result_1, 4)
  
  # Test with fac_rr = 2 scalers (0.98, 1.02)
  result_2 <- phytoclass:::Replace_Rand(Fmat_list, i, S, cm,
                                        min.scaler = 0.98, max.scaler = 1.02)
  expect_type(result_2, "list")
  expect_length(result_2, 4)
  
  # Test with fac_rr = 3 scalers (0.97, 1.03)
  result_3 <- phytoclass:::Replace_Rand(Fmat_list, i, S, cm,
                                        min.scaler = 0.97, max.scaler = 1.03)
  expect_type(result_3, "list")
  expect_length(result_3, 4)
})
