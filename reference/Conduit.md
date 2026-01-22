# Conduit between minimise_elements function and Fac_F_R of steepest descent algorithm.

Conduit between minimise_elements function and Fac_F_R of steepest
descent algorithm.

## Usage

``` r
Conduit(Fmat, place, S, cm, c_num = c(1, 2, 3))
```

## Arguments

- Fmat:

  A matrix of contribution values for each pigment and taxa pair

- place:

  Indices of elements to be modified

- S:

  Matrix of pigment sample measurements

- cm:

  Vector of weights for each column

- c_num:

  A numeric value (1, 2, or 3) to select which scaler values to use

## Value

A list containing three elements: `1`: The modified F matrix `2`: Number
of iterations or modifications made `3`: The original F matrix before
modifications

## Examples

``` r
Fmat <- as.matrix(phytoclass::Fm)
S <- as.matrix(phytoclass::Sm)
S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)

# Get F_initial from NNLS_MF as done in Steepest_Descent
F_initial <- phytoclass::NNLS_MF(Fmat, S, S_weights)

result <- phytoclass:::Conduit(F_initial[[1]], place, S, S_weights, c_num = 3)
```
