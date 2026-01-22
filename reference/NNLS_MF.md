# Performs the non-negative matrix factorisation for given phytoplankton pigments and pigment ratios, to attain an estimate of phytoplankton class abundances.

Performs the non-negative matrix factorisation for given phytoplankton
pigments and pigment ratios, to attain an estimate of phytoplankton
class abundances.

## Usage

``` r
NNLS_MF(Fn, S, S_weights = NULL)
```

## Arguments

- Fn:

  Pigment to Chl a matrix

- S:

  Sample data matrix – a matrix of pigment samples

- S_weights:

  Weights for each column

## Value

A list containing

1.  The F matrix (pigment: Chl *a*) ratios

2.  The root mean square error (RMSE)

3.  The C matrix (class abundances for each group)

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
