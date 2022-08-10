#' wrangles the data for the simulated annealing algorithm (maxs, mins etc...)
#'
#' @param Fl 
#'
#' @return
#' @export
#'
#' @examples
Wrangling <- function(Fl){
  Fd <- Fl
  Fmin <- as.matrix(Fd)   #### Set up Fmin matrix
  Fmin[Fmin>0] <- min  # set all non-zero elements to the minimum values (imported from csv )
  Fmax <- as.matrix(Fd)
  Fmax[Fmax>0] <- max
  chlv <- Fd[,ncol(Fd)]  ##### The chlorophyll values once weighted to rowsums for initial F matrix
  Fmin <- Fmin * chlv #### multiply the minimum value by weighted chlorophyll to updated ratios
  Fmin[,ncol(Fmin)] <- chlv  #### Reassign correct initial chl values
  Fmax <- Fmax * chlv  #### same for max values
  Fmax[,ncol(Fmax)] <- chlv
  
  Fmin <- vectorise(Fmin[,1:ncol(Fmin)-1]) #### vectorise function outputs all non-zero elements as a vector (excluding chl column)
  Fmax <- vectorise(Fmax[,1:ncol(Fmax)-1]) #### vectorise function outputs all non-zero elements as a vector (excluding chl column)
  SE <- vectorise(Fd[,1:ncol(Fd)-1]) #### vectorise function outputs all non-zero elements as a vector (excluding chl column)
  
  res <- list(Fmin, Fmax, SE,chlv)
  return(res)
}


