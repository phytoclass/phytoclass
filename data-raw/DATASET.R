## code to prepare `DATASET` dataset goes here

Sp <- read.csv("data-raw/Sp.csv", row.names = 1)
use_data(Sp, overwrite = TRUE)

Fp <- read.csv("data-raw/Fp.csv", row.names = 1)
use_data(Fp, overwrite = TRUE)


