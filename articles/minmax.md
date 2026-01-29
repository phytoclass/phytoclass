# Table of Default F MinMax Values

Below is the default table used to define the minimum and maximum values
for the pigment-to-chlorophyll-a ratios of each pigment-taxa pair. These
values are used to constrain the gradient decent algorithm.

``` r
print(phytoclass::min_max)
```

    ##                Class Pig_Abbrev    min   max
    ## 1                Syn        Zea 0.0800 2.500
    ## 2       Chlorophytes       Neox 0.0300 0.230
    ## 3       Chlorophytes       Viol 0.0110 0.470
    ## 4       Chlorophytes        Lut 0.0900 0.280
    ## 5       Chlorophytes        Zea 0.0020 0.099
    ## 6       Chlorophytes      Chl_b 0.1800 0.930
    ## 7      Prasinophytes       Neox 0.0180 0.150
    ## 8      Prasinophytes        Pra 0.0600 0.880
    ## 9      Prasinophytes       Viol 0.0080 0.150
    ## 10     Prasinophytes        Lut 0.0010 0.180
    ## 11     Prasinophytes        Zea 0.0150 0.350
    ## 12     Prasinophytes      Chl_b 0.4400 1.030
    ## 13      Cryptophytes       Allo 0.1600 0.790
    ## 14         Diatoms-1     Chl_c1 0.0030 0.057
    ## 15         Diatoms-1       Fuco 0.3000 1.700
    ## 16         Diatoms-2     Chl_c3 0.0460 0.270
    ## 17         Diatoms-2       Fuco 0.3300 1.230
    ## 18      Pelagophytes     Chl_c3 0.1500 0.310
    ## 19      Pelagophytes     X19but 0.2400 3.100
    ## 20      Pelagophytes       Fuco 0.3100 1.370
    ## 21 Dinoflagellates-1        Per 0.2900 1.600
    ## 22     Haptophytes-1     Chl_c1 0.0200 0.072
    ## 23     Haptophytes-1       Fuco 0.2000 0.340
    ## 24     Haptophytes-6     Chl_c3 0.0500 0.300
    ## 25     Haptophytes-6     X19but 0.0010 0.190
    ## 26     Haptophytes-6 ChlcMGDG18 0.0090 0.300
    ## 27     Haptophytes-6 ChlcMGDG14 0.0090 0.300
    ## 28     Haptophytes-6     X19hex 0.0300 1.500
    ## 29     Haptophytes-8     Chl_c3 0.0210 0.550
    ## 30     Haptophytes-8     X19but 0.0001 0.980
    ## 31     Haptophytes-8       Fuco 0.0110 2.150
    ## 32     Haptophytes-8     X19hex 0.0001 1.400
    ## 33       Haptophytes     Chl_c3 0.0210 0.550
    ## 34       Haptophytes     X19but 0.0001 0.980
    ## 35       Haptophytes       Fuco 0.0110 2.150
    ## 36       Haptophytes     X19hex 0.0001 1.400
    ## 37       Haptophytes ChlcMGDG18 0.0090 0.300
    ## 38       Haptophytes ChlcMGDG14 0.0090 0.300
    ## 39               Pro        Zea 0.0500 1.200
    ## 40               Syn      Tchla 1.0000 1.000
    ## 41      Chlorophytes      Tchla 1.0000 1.000
    ## 42     Prasinophytes      Tchla 1.0000 1.000
    ## 43      Cryptophytes      Tchla 1.0000 1.000
    ## 44         Diatoms-1      Tchla 1.0000 1.000
    ## 45         Diatoms-2      Tchla 1.0000 1.000
    ## 46      Pelagophytes      Tchla 1.0000 1.000
    ## 47 Dinoflagellates-1      Tchla 1.0000 1.000
    ## 48 Dinoflagellates-2      Tchla 1.0000 1.000
    ## 49     Haptophytes-1      Tchla 1.0000 1.000
    ## 50     Haptophytes-6      Tchla 1.0000 1.000
    ## 51       Haptophytes      Tchla 1.0000 1.000
    ## 52     Haptophytes-8      Tchla 1.0000 1.000
    ## 53               Pro     Dvchla 1.0000 1.000

These values are derived from the following paper among other sources:

> Mackey, M. D., Mackey, D. J., Higgins, H. W., & Wright, S. W. (1996).
> CHEMTAX—a program for estimating class abundances from chemical
> markers: application to HPLC measurements of phytoplankton. Marine
> Ecology Progress Series, 144(1/3), 265–283.
> <http://www.jstor.org/stable/24857270>
