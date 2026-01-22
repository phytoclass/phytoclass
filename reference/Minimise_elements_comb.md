# Part of the steepest descent algorithm

Part of the steepest descent algorithm

## Usage

``` r
Minimise_elements_comb(Fmat, place, S, cm, c1_num = c(1, 2, 3))
```

## Arguments

- Fmat:

  The F matrix

- place:

  A vector of indices indicating which elements to adjust

- S:

  A matrix of samples (rows) and pigments (columns)

- cm:

  A vector of bounded weights for each pigment

- c1_num:

  A numeric vector (1, 2, or 3) indicating which scaler values to use

## Value

A list containing the optimized F matrix with reduced error and vector
of indices where improvements were found

## Examples

``` r
 Fmat <- as.matrix(phytoclass::Fm)
 S <- as.matrix(phytoclass::Sm)
 S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
 place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)

 # Get F_initial from NNLS_MF as done in Steepest_Descent
 F_initial <- phytoclass::NNLS_MF(Fmat, S, S_weights)

 # Run Minimise_elements_comb with c1_num = 3 (as in Steepest_Descent)
 result <- phytoclass:::Minimise_elements_comb(
   F_initial[[1]], place, S, S_weights, c1_num = 3
 )
```
