# Normalize F matrix specifically for Prochlorococcus pigments

Normalizes pigment ratios differently for Prochlorococcus vs other
groups, using divinyl chlorophyll a for Prochlorococcus and chlorophyll
a for others.

## Usage

``` r
Prochloro_Normalise_F(Fmat)
```

## Arguments

- Fmat:

  Matrix with pigment ratios, where the last column is chlorophyll a and
  second-to-last column is divinyl chlorophyll a

## Value

A list containing: `1`: Normalized matrix where each row is scaled by
its biomass marker `2`: Vector of row sums from the scaled matrix

## Examples

``` r
# Create sample F matrix with Prochlorococcus
Fmat <- as.matrix(phytoclass::Fp)
result <- phytoclass:::Prochloro_Normalise_F(Fmat)
```
