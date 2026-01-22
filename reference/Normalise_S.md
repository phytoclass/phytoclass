# Normalise matrix to row sum

This function normalises each column in S to row sum

## Usage

``` r
Normalise_S(S)
```

## Arguments

- S:

  A matrix or data.frame to be normalized

## Value

A matrix

## Examples

``` r
# Create a sample matrix
S <- as.matrix(phytoclass::Sm)
normalized <- phytoclass:::Normalise_S(S)
```
