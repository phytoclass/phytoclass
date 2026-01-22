# Wrangle data to vectors

Converts data-types and selects data for randomisation in the simulated
annealing algorithm

## Usage

``` r
Wrangling(Fl, min.val, max.val)
```

## Arguments

- Fl:

  A matrix of the initial F matrix (i.e. pigment ratio matrix)

- min.val:

  A vector of the minimum values for each non-zero pigment ratios

- max.val:

  A vector of the maximum values for each non-zero pigment ratios

## Value

    A list containing following components:
    - A vector Fmin with the minimum pigment ratio values
    - A vector Fmax with the maximum pigment ratio values
    - A vector SE with the current pigment ratio values
    - A vector chlv with the pigment ratio values for the last column in Fl
