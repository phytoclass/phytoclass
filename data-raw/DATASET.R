## code to prepare `DATASET` dataset goes here

min_max <- read.csv("data-raw/Simon_pigments_new.csv")
use_data(min_max, overwrite = TRUE)

Sp <- read.csv("data-raw/Sp.csv", row.names = 1)
use_data(Sp, overwrite = TRUE)

Fp <- read.csv("data-raw/Fp.csv", row.names = 1)
use_data(Fp, overwrite = TRUE)


