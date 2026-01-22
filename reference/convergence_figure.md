# Convergence Figure

A figure to show the pigment ratios for each phytoplankton group for
each iteration.

## Usage

``` r
convergence_figure(fm_iter, niter = NULL)
```

## Arguments

- fm_iter:

  A data.frame with columns of iter, phyto, pigment and ratio

- niter:

  Optional: the number of iterations on the x axis. If `NULL`, will
  extract from the `iter` column of `fm_iter`.

## Value

A figure with each pigment ratio per iteration per group

## Examples

``` r
# ADD_EXAMPLES_HERE
```
