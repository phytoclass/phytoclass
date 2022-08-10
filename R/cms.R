#Ok, so these functions below are to add weights to the data.
# I want this to be an option to the user, with the default being no weights added. 

cms<-function(S){
  n <- colMeans(S)
  S <- n^-1
  S <- ifelse(S>50,50,S)
  S[length(S)] =1
  return(S)
}