# Example : Minimal

This vignette outlines the simplest possible usage of phytoclass’s
simulated annealing.

## Preparing the S Matrix

To run phytoclass you must prepare a sample matrix (aka the `S matrix`).
For this we recommend creating a comma-separated-value (`.csv`) file.
Several example files are included in the [`./vignettes/` directory of
the phytoclass source
code](https://github.com/phytoclass/phytoclass/tree/main/vignettes).
`.csv` files can be created in any plain text editor. Below is an
example S-matrix `.csv` file content:

                Per  X19but    Fuco    Neox     Pra    Viol  X19hex    Allo     Zea Lut ChlcMGDG18 ChlcMGDG14   Chl_b   Tchla
    Sample_1 0.0000 0.03024 0.06225 0.00557 0.01407 0.00759 0.08224 0.00188 0.00201   0          0          0 0.08661 0.45851
    Sample_2 0.0000 0.01084 0.02864 0.00111 0.00351 0.00144 0.01497 0.00144 0.00191   0          0          0 0.01473 0.14571
    Sample_3 0.0000 0.01560 0.21720 0.00640 0.00920 0.00000 0.01740 0.00580 0.00360   0          0          0 0.05700 0.61270
    Sample_4 0.0000 0.01770 0.23470 0.00700 0.01150 0.00000 0.01890 0.00540 0.00380   0          0          0 0.06190 0.62070
    Sample_5 0.0000 0.02520 0.29520 0.00990 0.01300 0.00000 0.02110 0.00140 0.00760   0          0          0 0.05780 0.53020
    Sample_6 0.0102 0.02220 0.22750 0.00760 0.01070 0.00000 0.01900 0.00000 0.00300   0          0          0 0.04530 0.40570

This example includes only 14 samples (rows), with the pigment
concentration of each sample included on each row. 14 is the minumum
number of rows required to apply simulated annealing. The first row
includes the names for each pigment. These **pigment abbreviations must
match exactly** to the pigments in the F-matrix. For this example we are
using the default `phytoclass:Fm` F matrix, and the column names in the
example above should be used.

> NOTE: Instead of using a plain text editor, `.csv` files can also be
> created using a spreadsheet editor (MS Excel, Google Docs, LibreOffice
> Calc, etc.) by using “save as” or “export” to `.csv`.

## Simulated Annealing

To run the simulated annealing call the function with your S matrix
passed as the only argument.

``` r
# Load the sample matrix from a file
# If your file does not contain sample names (e.g., Sample_1, Sample_2), you may omit `row.names = 1`
# In that case, do: S_matrix <- read.csv("vignettes/custom-example-S.csv")
S_matrix <- read.csv("vignettes/custom-example-S.csv", row.names = 1)

# run simulated annealing 
results <- phytoclass::simulated_annealing(S_matrix)
```

## Viewing Results

Results of the simulated annealing run are in the `results` object
returned by the function.

``` r
results$`condition number`
results$RMSE
results$MAE
results$Error
results$`Class abundances`
results$Figure
```
