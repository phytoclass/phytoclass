# Part of the steepest descent algorithm that works to reduce error given the S and F matrices.

Part of the steepest descent algorithm that works to reduce error given
the S and F matrices.

## Usage

``` r
Fac_F_RR(Fmat, vary, S, cm, fac_rr = c(1, 2, 3), place = NULL)
```

## Arguments

- Fmat:

  A list containing `F matrix`, `RMSE` and `C matrix`

- vary:

  Indices of non-zero elements to vary in the optimization

- S:

  A matrix of samples (rows) and pigments (columns)

- cm:

  A vector of bounded weights for each pigment

- fac_rr:

  A numeric value (1, 2, or 3) to select which scaler values to use: 1:
  (0.99, 1.01), 2: (0.98, 1.02), 3: (0.97, 1.03)

- place:

  A vector of all the indices of non-zero pigment ratios

## Value

A list containing two elements: `1`: Updated F matrix after optimization
`2`: Vector of indices

## Examples

``` r
 # Setup based on Minimise_elements_comb usage
 Fmat <- as.matrix(phytoclass::Fm)
 S <- as.matrix(phytoclass::Sm)
 cm <- as.numeric(phytoclass:::Bounded_weights(S))
 place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)

 # Get F.new from Conduit as done in Minimise_elements_comb
 f <- phytoclass:::Conduit(Fmat, place, S, cm, c_num = 3)
 F.new <- f[[1]]

 # Set place1 as done in Minimise_elements_comb
 place1 <- place # place1 = place when c1_num != 1

 # Run Fac_F_RR
 result <- phytoclass:::Fac_F_RR(F.new, vary = place, place = place1, S, cm, fac_rr = 3)
```
