# This is the main phytoclass algorithm. It performs simulated annealing algorithm for S and F matrices. See the examples (Fm, Sm) for how to set up matrices, and the vignette for more detailed instructions. Different pigments and phytoplankton groups may be used.

This is the main phytoclass algorithm. It performs simulated annealing
algorithm for S and F matrices. See the examples (Fm, Sm) for how to set
up matrices, and the vignette for more detailed instructions. Different
pigments and phytoplankton groups may be used.

## Usage

``` r
simulated_annealing(
  S,
  Fmat = NULL,
  user_defined_min_max = NULL,
  do_matrix_checks = TRUE,
  niter = 500,
  step = 0.009,
  weight.upper.bound = 30,
  verbose = TRUE,
  seed = NULL,
  check_converge = 100,
  alt_pro_name = NULL
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

  Set number to reproduce the same results

- check_converge:

  TRUE/FALSE/integer; set the number of F matrices to for convergence
  checking

- alt_pro_name:

  Optional: additional alternate versions of divinyl-chlorophyll-a
  spellings used to detect prochlorococcus (Default: "dvchl", "dvchla",
  "dv_chla")

## Value

A list containing

1.  Fmat matrix

2.  RMSE (Root Mean Square Error)

3.  condition number

4.  Class abundances

5.  Figure (plot of results)

6.  MAE (Mean Absolute Error)

7.  Error

8.  F_mat_iter

9.  converge_plot

## Examples

``` r
# Using the built-in matrices Sm and Fm
set.seed(5326)
sa.example <- simulated_annealing(Sm, Fm, niter = 5)
#> 
#> Condition number = 826
#> 
#> Iterations:         001 of 5 
#> Current error:      0.0269 
#> Neighbour's error:  0.0269 
#> Temperature (%):    99.1 
#> Iterations:         002 of 5 
#> Current error:      0.0255 
#> Neighbour's error:  0.0255 
#> Temperature (%):    98.21 
#> Iterations:         003 of 5 
#> Current error:      0.0255 
#> Neighbour's error:  0.0273 
#> Temperature (%):    97.32 
#> Iterations:         004 of 5 
#> Current error:      0.0255 
#> Neighbour's error:  0.027 
#> Temperature (%):    96.45 
#> Iterations:         005 of 5 
#> Current error:      0.0237 
#> Neighbour's error:  0.0237 
#> Temperature (%):    95.58 
sa.example$Figure
```
