# Randomise value by applying scaling factors within specified bounds. Small values (\< 0.001) are set to 0.001 to avoid numerical issues.

Randomise value by applying scaling factors within specified bounds.
Small values (\< 0.001) are set to 0.001 to avoid numerical issues.

## Usage

``` r
Randomise_elements(x, min.scaler, max.scaler)
```

## Arguments

- x:

  The element value to randomize

- min.scaler:

  Minimum scaling factor to apply (e.g., 0.99 for 1% decrease)

- max.scaler:

  Maximum scaling factor to apply (e.g., 1.01 for 1% increase)

## Value

A numeric value between x*min.scaler and x*max.scaler, rounded to 4
decimals

## Examples

``` r
# Randomize a single value
x <- 0.5
new_value <- phytoclass:::Randomise_elements(x, 0.99, 1.01)  # +/- 1% change

# Handle small values
small_x <- 0.0005
new_small <- phytoclass:::Randomise_elements(small_x, 0.99, 1.01)  # will use 0.001
```
