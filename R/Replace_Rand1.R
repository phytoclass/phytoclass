# Apply randomization functions to matrices
Replace_Rand1 <- function(F,i){
  F.new <- as.matrix(replace(F[[1]],i,rand1(F[[1]][i]))) # randomise first element of matrix
  F.new <- Fac_F(F.new)
  v <- which(F.new[[2]] < F[[2]]) #Which elements decrease the error? Store the location of the elements that decrease it 
  res <- c(F.new,v) 
  return(res)
}
