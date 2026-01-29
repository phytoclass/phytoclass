# Changelog

## phytoclass 2.3.1

### Enhancements

- Updated `min_max.rda` using new pigment ratio min/max from Simon
  Wright

## phytoclass 2.3.0

This version update removes internal functions that became unused once
combined into a single function. The sets of functions include
“Conduit”, “Fac_F_RR”, “Minimise_elements” and “Random_neighbor”. \##
New features \*
[`convergence_figure()`](../reference/convergence_figure.md) was added
in
[`simulated_annealing_Prochloro()`](../reference/simulated_annealing_Prochloro.md)
\* Progress bar was added in
[`simulated_annealing_Prochloro()`](../reference/simulated_annealing_Prochloro.md)
\* [`Cluster()`](../reference/Cluster.md) now accepts a `row_ids` input
vector to name each sample \* [`Cluster()`](../reference/Cluster.md)
exports an `assigments` item that relates the ID to the cluster \*
[`Cluster()`](../reference/Cluster.md) checks if data.frame input has a
column of strings and removes them prior to clustering

### Bug Fixes

- The figure functions would `error` if the number of colors needed
  exceeded 11. Now, the figure functions can handle additional colors.
- [`Cluster()`](../reference/Cluster.md) had an issue if given a matrix
  without rownames and would error, now will convert to data.frame prior
  to running boxcar normalising

## phytoclass 2.2.0

- Add automatic switch to Prochlorococcus function based on dvchl column
- added convergence plot

## phytoclass 2.1.0

- performance improvements
- improved documentation
- - continuous integration
- - tests
- improved matrix checks
- improved error messages

## phytoclass 2.0.0

CRAN release: 2024-11-14

## phytoclass 1.2.0

CRAN release: 2024-05-11

## phytoclass 1.1.0

## phytoclass 0.0.0.9000

- Added a `NEWS.md` file to track changes to the package.

## phytoclass 1.1.0

### Breaking changes

### New features

We have added a function which allows for chlorophyll derivations for
prochlorococcus.

### Minor improvements and fixes

We have improved documentation, and relabelled groups, such as
‘Diatom-1’ instead of ‘Diatom-A’.
