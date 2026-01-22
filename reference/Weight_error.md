# Apply weights to F/S matrices by diagonal multiplication

Apply weights to F/S matrices by diagonal multiplication

## Usage

``` r
Weight_error(S, cm)
```

## Arguments

- S:

  Matrix to be weighted

- cm:

  Vector of weights to be applied to columns of S

## Value

A matrix with weighted columns (S %\*% diag(cm))

## Examples

``` r
# Create sample matrix and weights
S <- as.matrix(phytoclass::Sm)
Fmat <- as.matrix(phytoclass::Fm)
S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
weighted <- phytoclass:::Weight_error(Fmat, S_weights)
```
