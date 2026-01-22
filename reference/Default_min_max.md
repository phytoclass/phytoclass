# Sets the default minimum and maximum values for phytoplankton groups pigment ratios. To use this function, pigment and phytoplankton group names will need to fit the naming criteria of phytoclass.

Sets the default minimum and maximum values for phytoplankton groups
pigment ratios. To use this function, pigment and phytoplankton group
names will need to fit the naming criteria of phytoclass.

## Usage

``` r
Default_min_max(min_max, Fmat)
```

## Arguments

- min_max:

  A data.frame with 4 columns for class, pigment, min and max values

- Fmat:

  F matrix with phytoplankton groups as rows and pigments as columns

## Value

A list containing two elements: `1`: Vector of minimum values for each
non-zero pigment ratio `2`: Vector of maximum values for each non-zero
pigment ratio

## Examples

``` r
# Create a sample F matrix
Fmat <- phytoclass::Fm

# Create min_max data frame
min_max <- phytoclass::min_max
result <- phytoclass:::Default_min_max(min_max, Fmat)
```
