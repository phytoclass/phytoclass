# Add weights to the data, bound at a maximum.

Add weights to the data, bound at a maximum.

## Usage

``` r
Bounded_weights(S, weight.upper.bound = 30)
```

## Arguments

- S:

  Sample data matrix – a matrix of pigment samples

- weight.upper.bound:

  Upper bound for weights (default is 30)

## Value

A vector with upper bounds for weights

## Examples

``` r
Bounded_weights(Sm, weight.upper.bound = 30)
#>        Per     X19but       Fuco       Neox        Pra       Viol     X19hex 
#>  30.000000  30.000000   3.405156  30.000000  30.000000  30.000000  15.742517 
#>       Allo        Zea        Lut ChlcMGDG18 ChlcMGDG14      Chl_b      Tchla 
#>  30.000000  30.000000  30.000000  30.000000  30.000000   9.185585   1.000000 
```
