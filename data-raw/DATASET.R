## code to prepare `DATASET` dataset goes here

# row names are in the first column
S <- read.csv("data-raw/S_example.csv", header = TRUE, row.names = 1)
# remove column at end containing "Shelf"
Sm <- S[, -16]
usethis::use_data(Sm, overwrite = TRUE)

Fm <- read.csv("data-raw/F_example.csv", header = TRUE, row.names = 1)
usethis::use_data(Fm, overwrite = TRUE)

min_max <-read.csv("data-raw/Circumpolar_minmax.csv", header = TRUE) # Read pig:chl inital
min_max <- min_max[, c("Class", "Pig_Abbrev", "min", "max")]
usethis::use_data(min_max , overwrite = TRUE)
