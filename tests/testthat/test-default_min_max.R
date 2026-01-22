test_that("Default_min_max returns correct min and max for numeric vector", {
  # Example from the documentation
  Fmat <- phytoclass::Fm
  min_max <- phytoclass::min_max
  result <- phytoclass:::Default_min_max(min_max, Fmat)
  
  # Should return a list with min and max vectors
  expect_type(result, "list")
  expect_length(result, 2)
  expect_true(is.numeric(result[[1]]))
  expect_true(is.numeric(result[[2]]))
  expect_equal(length(result[[1]]), length(result[[2]]))
})
