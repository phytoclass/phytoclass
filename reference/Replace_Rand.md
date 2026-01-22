# Select the new F matrix element with lowest error in the steepest descent algorithm. Randomly modifies a single element and checks if the modification reduces the error.

Select the new F matrix element with lowest error in the steepest
descent algorithm. Randomly modifies a single element and checks if the
modification reduces the error.

## Usage

``` r
Replace_Rand(Fmat, i, S, cm, min.scaler, max.scaler)
```

## Arguments

- Fmat:

  A list containing the F matrix, RMSE, and other components

- i:

  The index of the element to modify in the F matrix

- S:

  Sample data matrix - matrix of pigment samples

- cm:

  A vector of bounded weights for each pigment

- min.scaler:

  Minimum scaling factor to apply (e.g., 0.99 for 1% decrease)

- max.scaler:

  Maximum scaling factor to apply (e.g., 1.01 for 1% increase)

## Value

A list containing:

- F matrix:

  The modified F matrix

- RMSE:

  Root mean square error of the new solution

- C matrix:

  The concentration matrix

- Improved:

  Logical indicating if the modification reduced error

## Examples

``` r
 # Setup based on Fac_F_RR usage
 Fmat <- as.matrix(phytoclass::Fm)
 S <- as.matrix(phytoclass::Sm)
 cm <- as.numeric(phytoclass:::Bounded_weights(S))

 # Get Fmat as a list from NNLS_MF (as used in Fac_F_RR)
 Fmat_list <- phytoclass::NNLS_MF(Fmat, S, cm)

 # Test with a single index
 i <- 1 # first non-zero element to modify
 min.scaler <- 0.99
 max.scaler <- 1.01

 # Run Replace_Rand
 result <- phytoclass:::Replace_Rand(Fmat_list, i, S, cm, min.scaler, max.scaler)
```
