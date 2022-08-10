
# This returns all elements that aren't 0 in the F matrix
vectorise <- function(Fg){
  g <- vector()
  for (i in Fg){
    if (i > -0){
      g[[length(g)+1]] <- i
      
    }
    
  }
  return(g)
}
