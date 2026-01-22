# Perform matrix factorisation for phytoplankton pigments and pigments ratios

Performs the non-negative matrix factorisation for given phytoplankton
pigments and pigment ratios, to attain an estimate of phytoplankton
class abundances.

## Usage

``` r
NNLS_MF_Final(Fn, S, S_Chl, S_weights, S_dvChl = NULL)
```

## Arguments

- Fn:

  F matrix with pigment ratios for each phytoplankton class

- S:

  Sample data matrix of pigment measurements

- S_Chl:

  Vector of chlorophyll a concentrations for each sample

- S_weights:

  Vector of weights for each pigment

- S_dvChl:

  Optional vector of divinyl chlorophyll concentrations for
  Prochlorococcus

## Value

A list containing the following elements:

- F matrix:

  The normalized F matrix of pigment ratios

- RMSE:

  Root mean square error of the fit

- condition number:

  Condition number of Fn %\*% t(S)

- Class abundances:

  Data frame of phytoplankton class abundances

- Figure:

  Plot of the results

- MAE:

  Mean absolute error for each pigment

- Error:

  Residual error matrix

## Details

Unlike NNLS_ML(), it also removes any weighting and normalisation, and
also multiplies relative abundances by chlorophyll values to determine
the biomass of phytoplankton groups.

## Examples

``` r
 Fmat <- as.matrix(phytoclass::Fm)
 S <- as.matrix(phytoclass::Sm)
 S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
 S_Chl <- S[, ncol(S)]
 # Run NNLS_MF_Final
 result <- phytoclass:::NNLS_MF_Final(Fmat, S, S_Chl, S_weights)
```
