# Perform simulated annealing algorithm for samples with divinyl chlorophyll and prochlorococcus. Chlorophyll must be the final column of both S and F matrices, with Divinyl Chlorophyll a the 2nd to last column. See how the example Sp and Fp matrices are organised.

Perform simulated annealing algorithm for samples with divinyl
chlorophyll and prochlorococcus. Chlorophyll must be the final column of
both S and F matrices, with Divinyl Chlorophyll a the 2nd to last
column. See how the example Sp and Fp matrices are organised.

## Usage

``` r
simulated_annealing_Prochloro(
  S,
  Fmat = NULL,
  user_defined_min_max = NULL,
  do_matrix_checks = TRUE,
  niter = 500,
  step = 0.009,
  weight.upper.bound = 30,
  verbose = TRUE,
  seed = NULL,
  check_converge = 100
)
```

## Arguments

- S:

  Sample data matrix – a matrix of pigment samples

- Fmat:

  Pigment to Chl a matrix

- user_defined_min_max:

  data frame with some format as min_max built-in data

- do_matrix_checks:

  This should only be set to TRUE when using the default values. This
  will remove pigment columns that have column sums of 0. Set to FALSE
  if using customised names for pigments and phytoplankton groups

- niter:

  Number of iterations (default is 500)

- step:

  Step ratio used (default is 0.009)

- weight.upper.bound:

  Upper limit of the weights applied (default value is 30).

- verbose:

  Logical value. Output error and temperature at each iteration. Default
  value of TRUE

- seed:

  Set seed number to reproduce the same results

- check_converge:

  TRUE/FALSE/integer; set the number of F matrices to for convergence
  checking

## Value

A list containing

1.  Fmat matrix

2.  RMSE (Root Mean Square Error)

3.  condition number

4.  Class abundances

5.  Figure (plot of results)

6.  MAE (Mean Absolute Error)

7.  Error

## Examples

``` r
# Using the built-in matrices Sp and Fp.
set.seed(5326)
sa.example <- simulated_annealing_Prochloro(Sp, Fp, niter = 1)
#> 
#> Condition number = 11601
#> 
#> Iterations:         001 of 1 
#> Current error:      0.009 
#> Neighbour's error:  0.009 
#> Temperature (%):    99.1 
#> Using 15 color palette for `convergence plot`.
sa.example$Figure
```
