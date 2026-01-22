# Common Issues : F Matrix

## non-conformable arguments

> Error in Fn %\*% t(S) : non-conformable arguments

This error can mean there is some disagreement between the S and F
matricies.

## non-numeric argument to binary operator

> non-numeric argument to binary operator

This error can mean that the F matrix has extra columns that should not
be there. \* [ref](https://github.com/phytoclass/phytoclass/issues/10)

## x needs two+ dimensions

> ‘x’ must be an array of at least two dimension

<https://github.com/phytoclass/phytoclass/issues/11>

## subscript out of bounds:

> Error in x\[subset & !is.na(subset), vars, drop = drop\] : subscript
> out of bounds

<https://github.com/phytoclass/phytoclass/issues/12>

This error can happen when the sample names are included in the S
matrix. Current version of phytoclass does not support sample names.

## ‘x’ must be numeric

> Error in rowSums(Fmat) : ‘x’ must be numeric

This can mean that the F matrix has a numeric index (1,2,3,etc) column
that must be removed. This is sometimes introduced by `read.csv` or
`readRDS`. A workaround for this:

``` r
# === remove numeric rownames introduced by read.csv
if (all(grepl("^[0-9]+$", rownames(F_matrix)))) {
  print("dropping unneeded numeric index")
  # Set the first column as row names
  rownames(F_matrix) <- F_matrix[[1]]
  
  # Remove the first column
  F_matrix <- F_matrix[, -1] 
}
```
