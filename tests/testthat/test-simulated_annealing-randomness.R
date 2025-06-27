test_that("simulated_annealing output changes with different seeds", {
  Sm <- phytoclass::Sm
  Fm <- phytoclass::Fm
  
  out1 <- simulated_annealing(Sm, Fm, niter = 10, seed = 1234, verbose = FALSE)
  out2 <- simulated_annealing(Sm, Fm, niter = 10, seed = 4567, verbose = FALSE)
  
  # Compare class abundances
  expect_false(
    isTRUE(all.equal(out1$`Class abundances`, out2$`Class abundances`, tolerance = 1e-10))
  )
  # Compare RMSE values
  expect_false(
    isTRUE(all.equal(out1$RMSE, out2$RMSE, tolerance = 1e-10))
  )
  
  # Compare MAE vector
  expect_false(
    isTRUE(all.equal(out1$MAE, out2$MAE, tolerance = 1e-10))
  )
  
})
