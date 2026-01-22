# Performs the steepest descent algorithm for a set number of iterations to optimize the F matrix of pigment ratios.

Performs the steepest descent algorithm for a set number of iterations
to optimize the F matrix of pigment ratios.

## Usage

``` r
Steepest_Descent(Fmat, place, S, S_weights, num.loops)
```

## Arguments

- Fmat:

  Initial F matrix containing pigment ratios

- place:

  Vector of indices where F matrix has non-zero values

- S:

  Matrix of sample measurements (rows) and pigments (columns)

- S_weights:

  Vector of weights for each pigment in NNLS optimization

- num.loops:

  Maximum number of iterations to perform optimization

## Value

A list containing: `1`: The optimized F matrix `2`: Final RMSE value
`3`: The C matrix (class abundances for each group)

## Examples

``` r
 Fmat <- as.matrix(phytoclass::Fm)
 S <- as.matrix(phytoclass::Sm)
 S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
 place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
 num.loops <- 2
 # Run Steepest_Descent
 result <- phytoclass:::Steepest_Descent(Fmat, place, S, S_weights, num.loops)
```
