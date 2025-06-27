test_that("simulated_annealing output is the same with same seeds", {
  Sm <- phytoclass::Sm
  Fm <- phytoclass::Fm
  
  out1 <- simulated_annealing(Sm, Fm, niter = 10, seed = 1234, verbose = FALSE)
  out2 <- simulated_annealing(Sm, Fm, niter = 10, seed = 1234, verbose = FALSE)
  
  # These should match exactly
  expect_true(
    isTRUE(all.equal(out1$`Class abundances`, out2$`Class abundances`, tolerance = 1e-10))
  )
  expect_true(
    isTRUE(all.equal(out1$RMSE, out2$RMSE, tolerance = 1e-10))
  )
  expect_true(
    isTRUE(all.equal(out1$MAE, out2$MAE, tolerance = 1e-10))
  )
})
