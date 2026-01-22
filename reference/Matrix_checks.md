# Remove any column values that average 0. Further to this, also remove phytoplankton groups from the F matrix if their diagnostic pigment isn’t present.

Remove any column values that average 0. Further to this, also remove
phytoplankton groups from the F matrix if their diagnostic pigment isn’t
present.

## Usage

``` r
Matrix_checks(S, Fmat)
```

## Arguments

- S:

  Sample data matrix – a matrix of pigment samples

- Fmat:

  Pigment to Chl a matrix

## Value

Named list with new S and Fmat matrices

## Examples

``` r
MC <- Matrix_checks(Sm, Fm)  
Snew <- MC$Snew
```
