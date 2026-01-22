# Select a random neighbour when the previous random neighbour is beyond the minimum or maximum value. Part of the simulated annealing algorithm.

Select a random neighbour when the previous random neighbour is beyond
the minimum or maximum value. Part of the simulated annealing algorithm.

## Usage

``` r
Random_neighbour(f_new, Temp, chlv, N, place, S, S_weights, minF, maxF)
```

## Arguments

- f_new:

  Current F matrix of pigment ratios

- Temp:

  Current temperature in the annealing process

- chlv:

  Chlorophyll-a column to append to the matrix

- N:

  Indices of pigment ratios to be modified

- place:

  Vector of indices where values are non-zero in the F matrix

- S:

  Matrix of samples (rows) and pigments (columns)

- S_weights:

  Vector of weights for each pigment in NNLS

- minF:

  Minimum bounds for each pigment ratio

- maxF:

  Maximum bounds for each pigment ratio

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
 # Setup based on simulated_annealing usage
 Fmat <- as.matrix(phytoclass::Fm)
 S <- as.matrix(phytoclass::Sm)
 S_weights <- as.numeric(phytoclass:::Bounded_weights(S))
 place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
 min_max <- phytoclass::min_max
 minF <- min_max[[3]][seq_along(place)]
 maxF <- min_max[[4]][seq_along(place)]
 chlv <- rep(1, nrow(Fmat)) # typical usage in simulated_annealing
 Temp <- 0.5
 N <- place
 # Run Random_neighbour
 result <- phytoclass:::Random_neighbour(Fmat, Temp, chlv, N, place, S, S_weights, minF, maxF)
```
