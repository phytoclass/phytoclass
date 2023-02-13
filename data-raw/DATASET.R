## code to prepare `DATASET` dataset goes here

S <-read.csv("data-raw/Clusters_213.csv",header = TRUE,row.names = "X")
S <- S[ , 14:length(S)-1]
S <- S[, -c(18,11,2,1)]
S <- S[, !("V31" %in% names(S))]
usethis::use_data(S, overwrite = TRUE)

F <-read.csv("data-raw/Ratios_5.csv",header = TRUE,row.names = 'X') # Read pig:chl inital 
F <- F[,-1]
usethis::use_data(F, overwrite = TRUE)

min_max <-read.csv("data-raw/Circumpolar_minmax.csv", header = TRUE) # Read pig:chl inital 
usethis::use_data(min_max , overwrite = TRUE)
