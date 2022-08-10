Weight_error <- function(S){
  S <- S%*%diag(cm)
  return(S)
}