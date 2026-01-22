# Stand-alone version of steepest descent algorithm. This is similar to the CHEMTAX steepest descent algorithm. It is not required to use this function, and as results are not bound by minimum and maximum, results may be unrealistic.

Stand-alone version of steepest descent algorithm. This is similar to
the CHEMTAX steepest descent algorithm. It is not required to use this
function, and as results are not bound by minimum and maximum, results
may be unrealistic.

## Usage

``` r
Steepest_Desc(Fmat, S, num.loops)
```

## Arguments

- Fmat:

  Pigment to Chl a matrix

- S:

  Sample data matrix – a matrix of pigment samples

- num.loops:

  Number of loops/iterations to perform (no default)

## Value

A list containing

1.  The F matrix (pigment: Chl *a*) ratios

2.  RMSE (Root Mean Square Error)

3.  Condition number

4.  class abundances

5.  Figure (plot of results)

6.  MAE (Mean Absolute Error)

## Examples

``` r
MC <- Matrix_checks(Sm,Fm)
Snew <- MC$Snew
Fnew <- MC$Fnew
SDRes <- Steepest_Desc(Fnew,Snew, num.loops = 20)
```
