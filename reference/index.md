# Package index

## All functions

- [`Bounded_weights()`](Bounded_weights.md) : Add weights to the data,
  bound at a maximum.
- [`Cluster()`](Cluster.md) : Cluster things
- [`Fm`](Fm.md) : Fm data
- [`Fp`](Fp.md) : Fp data
- [`Matrix_checks()`](Matrix_checks.md) : Remove any column values that
  average 0. Further to this, also remove phytoplankton groups from the
  F matrix if their diagnostic pigment isn’t present.
- [`NNLS_MF()`](NNLS_MF.md) : Performs the non-negative matrix
  factorisation for given phytoplankton pigments and pigment ratios, to
  attain an estimate of phytoplankton class abundances.
- [`Sm`](Sm.md) : Sm data
- [`Sp`](Sp.md) : Sp data
- [`Steepest_Desc()`](Steepest_Desc.md) : Stand-alone version of
  steepest descent algorithm. This is similar to the CHEMTAX steepest
  descent algorithm. It is not required to use this function, and as
  results are not bound by minimum and maximum, results may be
  unrealistic.
- [`convergence_figure()`](convergence_figure.md) : Convergence Figure
- [`min_max`](min_max.md) : min_max data
- [`phyto_figure()`](phyto_figure.md) : Phytoplankton Class Abundance
  Figure
- [`simulated_annealing()`](simulated_annealing.md) : This is the main
  phytoclass algorithm. It performs simulated annealing algorithm for S
  and F matrices. See the examples (Fm, Sm) for how to set up matrices,
  and the vignette for more detailed instructions. Different pigments
  and phytoplankton groups may be used.
- [`simulated_annealing_Prochloro()`](simulated_annealing_Prochloro.md)
  : Perform simulated annealing algorithm for samples with divinyl
  chlorophyll and prochlorococcus. Chlorophyll must be the final column
  of both S and F matrices, with Divinyl Chlorophyll a the 2nd to last
  column. See how the example Sp and Fp matrices are organised.
