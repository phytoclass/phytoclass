# Table of Default F MinMax Values

Below is the default table used to define the minimum and maximum values
for the pigment-to-chlorophyll-a ratios of each pigment-taxa pair. These
values are used to constrain the gradient decent algorithm.

``` r
print(phytoclass::min_max)
```

    ##                Class Pig_Abbrev    min    max
    ## 1                Syn        Zea 0.0800 1.2123
    ## 2       Chlorophytes       Neox 0.0091 0.1085
    ## 3       Chlorophytes       Viol 0.0150 0.5875
    ## 4       Chlorophytes        Lut 0.1071 0.2295
    ## 5       Chlorophytes        Zea 0.0063 0.0722
    ## 6       Chlorophytes      Chl_b 0.1666 0.9254
    ## 7      Prasinophytes       Neox 0.0268 0.1026
    ## 8      Prasinophytes        Pra 0.0642 0.4369
    ## 9      Prasinophytes       Viol 0.0087 0.0996
    ## 10     Prasinophytes        Lut 0.0250 0.0669
    ## 11     Prasinophytes        Zea 0.0151 0.1396
    ## 12     Prasinophytes      Chl_b 0.4993 0.9072
    ## 13      Cryptophytes       Allo 0.2118 0.5479
    ## 14         Diatoms-1     Chl_c1 0.0021 0.0452
    ## 15         Diatoms-1       Fuco 0.3315 0.9332
    ## 16         Diatoms-2     Chl_c3 0.0189 0.1840
    ## 17         Diatoms-2       Fuco 0.3315 0.9332
    ## 18      Pelagophytes     Chl_c3 0.1471 0.2967
    ## 19      Pelagophytes     X19but 0.2457 1.4339
    ## 20      Pelagophytes       Fuco 0.3092 1.2366
    ## 21 Dinoflagellates-1        Per 0.3421 0.8650
    ## 22       Haptophytes     Chl_c3 0.0500 0.3500
    ## 23       Haptophytes     X19but 0.0819 0.2872
    ## 24       Haptophytes ChlcMGDG18 0.0090 0.3000
    ## 25       Haptophytes ChlcMGDG14 0.0090 0.3000
    ## 26       Haptophytes     X19hex 0.2107 1.3766
    ## 27               Pro        Zea 0.0800 1.2123
    ## 28       Haptophytes       Fuco 0.0090 0.4689
    ## 29       Haptophytes      Tchla 1.0000 1.0000
    ## 30         Diatoms-2      Tchla 1.0000 1.0000
    ## 31         Diatoms-1      Tchla 1.0000 1.0000
    ## 32      Cryptophytes      Tchla 1.0000 1.0000
    ## 33     Prasinophytes      Tchla 1.0000 1.0000
    ## 34      Chlorophytes      Tchla 1.0000 1.0000
    ## 35               Syn      Tchla 1.0000 1.0000
    ## 36 Dinoflagellates-1      Tchla 1.0000 1.0000
    ## 37 Dinoflagellates-2      Tchla 1.0000 1.0000
    ## 38      Pelagophytes      Tchla 1.0000 1.0000
    ## 39               Pro     Dvchla 1.0000 1.0000

These values are derived from the following paper among other sources:

> Mackey, M. D., Mackey, D. J., Higgins, H. W., & Wright, S. W. (1996).
> CHEMTAX—a program for estimating class abundances from chemical
> markers: application to HPLC measurements of phytoplankton. Marine
> Ecology Progress Series, 144(1/3), 265–283.
> <http://www.jstor.org/stable/24857270>
