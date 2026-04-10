# This function ensures S and F matrices are properly formatted and ordered for the simulated annealing function.

Some checks applied:

- drops columns with 0 values

- drops taxa with missing major pigments, which are indicated with a '2'

- drops pigments with \< 1% in samples

## Usage

``` r
Matrix_checks(S, Fmat)
```

## Arguments

- S:

  Sample data matrix – a matrix of pigment samples

- Fmat:

  Pigment to taxa matrix

## Value

Named list with new S and Fmat matrices

## Examples

``` r
MC <- Matrix_checks(Sm, Fm)  
Snew <- MC$Snew
```
