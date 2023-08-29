#' Converts data-types and selects data for randomisation in 
#' the simulated annealing algorithm
#' 
#' @keywords internal
#'
#' @param Fl   xx
#' @param min.val xx
#' @param max.val  xx
#'
#' @return
#'
#' @examples
Wrangling <- function(Fl, min.val, max.val){
  Fd <- Fl
  Fmin <- as.matrix(Fd)   #### Set up Fmin matrix
  Fmin <- Fmin[,1:ncol(Fmin)-1]
  Fmin[Fmin>0] <- min.val  # set all non-zero elements to the minimum values (imported from csv )
  Fmax <- as.matrix(Fd)
  Fmax <- Fmax[,1:ncol(Fmax)-1]
  Fmax[Fmax>0] <- max.val
  chlv <- Fd[,ncol(Fd)]  ##### The chlorophyll values once weighted to rowsums for initial F matrix
  Fmin <- Fmin * chlv #### multiply the minimum value by weighted chlorophyll to updated ratios
  Fmin <- cbind(Fmin,chlv)  #### Reassign correct initial chl values
  Fmax <- Fmax * chlv  #### same for max values
  Fmax <- cbind(Fmax,chlv)
  
  Fmin <- vectorise(Fmin[,1:ncol(Fmin)-1]) #### vectorise function outputs all non-zero elements as a vector (excluding chl column)
  Fmax <- vectorise(Fmax[,1:ncol(Fmax)-1]) #### vectorise function outputs all non-zero elements as a vector (excluding chl column)
  SE <- vectorise(Fd[,1:ncol(Fd)-1]) #### vectorise function outputs all non-zero elements as a vector (excluding chl column)
  
  res <- list(Fmin, Fmax, SE,chlv)
}

