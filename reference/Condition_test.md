# Calculate the mean condition number for randomized F matrices

Performs multiple simulations with randomized F matrices within given
bounds to assess the numerical stability of the system.

## Usage

``` r
Condition_test(S, Fn, min.val = NULL, max.val = NULL)
```

## Arguments

- S:

  Sample matrix of pigment measurements

- Fn:

  Initial F matrix of pigment ratios

- min.val:

  Optional vector of minimum values for each non-zero pigment ratio

- max.val:

  Optional vector of maximum values for each non-zero pigment ratio

## Value

Numeric value representing the mean condition number from 1000
simulations

## Examples

``` r
# Create sample matrices
Fmat <- as.matrix(phytoclass::Fm)
S <- as.matrix(phytoclass::Sm)
min_max <- phytoclass::min_max
min.val <- min_max[[3]]
max.val <- min_max[[4]]
# Use only non-chla columns for Fmat and S, as in simulated_annealing
Fmat_sub <- Fmat[, -ncol(Fmat)]
S_sub <- S[, -ncol(S)]
# Calculate mean condition number
cond <- phytoclass:::Condition_test(S_sub, Fmat_sub, min.val, max.val)
```
