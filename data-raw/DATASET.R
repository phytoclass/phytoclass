## code to prepare `DATASET` dataset goes here

# row names are in the first column
S <- read.csv("data-raw/S_example.csv", header = TRUE, row.names = 1)
# remove column at end containing "Shelf"
S <- S[, -16]
usethis::use_data(S, overwrite = TRUE)

F <- read.csv("data-raw/F_example.csv", header = TRUE, row.names = 1)
usethis::use_data(F, overwrite = TRUE)

min_max <-read.csv("data-raw/Circumpolar_minmax.csv", header = TRUE) # Read pig:chl inital 
usethis::use_data(min_max , overwrite = TRUE)
