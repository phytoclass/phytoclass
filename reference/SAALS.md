# Apply the steepest descent algorithm to optimize pigment ratios in phytoplankton classification. Loosely wraps Seepest_Descent function.

Apply the steepest descent algorithm to optimize pigment ratios in
phytoplankton classification. Loosely wraps Seepest_Descent function.

## Usage

``` r
SAALS(Ft, min.value, max.value, place, S, cm, num.loops)
```

## Arguments

- Ft:

  Initial F matrix containing pigment ratios

- min.value:

  UNUSED?

- max.value:

  UNUSED?

- place:

  Vector of indices where F matrix has non-zero values

- S:

  Matrix of sample measurements

- cm:

  Vector of bounded weights for each pigment

- num.loops:

  Maximum number of iterations for the steepest descent

## Value

A list containing: `1`: The optimized F matrix `2`: Final RMSE value

## Examples

``` r
 Fmat <- as.matrix(phytoclass::Fm)
 S <- as.matrix(phytoclass::Sm)
 S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
 place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
 num.loops <- 2
 # Run SAALS
 result <- phytoclass:::SAALS(Fmat, NULL, NULL, place, S, S_weights, num.loops)
```
