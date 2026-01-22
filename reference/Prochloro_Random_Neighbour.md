# Selects a random neighbour for a subset of non-zero pigments that are outside the min and max bounds for the simulated annealing algorithm, specifically handling Prochlorococcus pigments.

Selects a random neighbour for a subset of non-zero pigments that are
outside the min and max bounds for the simulated annealing algorithm,
specifically handling Prochlorococcus pigments.

## Usage

``` r
Prochloro_Random_Neighbour(
  f_new,
  Temp,
  chlv,
  chlvp,
  N,
  place,
  S,
  S_weights,
  minF,
  maxF
)
```

## Arguments

- f_new:

  F matrix of pigment ratios

- Temp:

  Temperature of the annealing

- chlv:

  Chlorophyll-a column

- chlvp:

  Dvchla column

- N:

  Indexs of pigment ratios to be changed

- place:

  Indexes in F matrix where values are non-zero

- S:

  S matrix of samples

- S_weights:

  Weights for NNLS algorithm

- minF:

  Minimum bounds for each phytoplankton group and pigments

- maxF:

  Maximum bounds for each phytoplankton group and pigments

## Value

A list containing:

- F matrix:

  The new F matrix with randomly modified values

- RMSE:

  Root mean square error of the new solution

- C matrix:

  The concentration matrix from NNLS

## Examples

``` r
 Fmat <- as.matrix(phytoclass::Fp)
 S <- as.matrix(phytoclass::Sp)
 S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
 place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)

 # Get min_max from package data
 min_max_mat <- phytoclass:::Default_min_max(phytoclass::min_max, Fmat[, -ncol(Fmat)])

 # Get bounds using Prochloro_Wrangling as in simulated_annealing_Prochloro
 f_c <- Fmat
 W0 <- phytoclass:::Prochloro_Wrangling(f_c, min_max_mat[[1]], min_max_mat[[2]])
 minF <- W0[[1]]
 maxF <- W0[[2]]

 # Extract chlv and chlvp from the matrix
 chlv <- f_c[, ncol(f_c)] # Chl a (Tchla)
 chlvp <- f_c[, ncol(f_c) - 1] # dvChl a

 Temp <- 0.5
 N <- place

 # Run Prochloro_Random_Neighbour
 result <- phytoclass:::Prochloro_Random_Neighbour(
   f_c, Temp, chlv, chlvp, N, place, S, S_weights, minF, maxF
 )
```
