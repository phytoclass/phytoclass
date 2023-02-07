#' Remove any column values that average 0. Further to this, also remove
#' phytoplankton groups from the F matrix if their diagnostic pigment
#' isn’t present. 
#'
#' @param S
#' @param F
#'
#' @return
#' @export
#'
#' @examples


Homogenise_matrices <- function(S, F){
  ba <- rownames(F)
  ba1<- which(ba =="Syn")
  if (unique(S$V31) =="Shelf" | unique(S$V31) =="SSIZ"){
    F <- F[-ba1,]
  }
  S <- S[,1:ncol(S)-1]
  b <- colSums(S)
  g <- mean(S$Tchla)
  ba <- rownames(F)
  ba1<- which(ba =="Syn")
  
  b <- nrow(S)
  c <- which(b ==0)
  g <- colSums(S != 0)
  l <- which(g/b <=.2)
  if(length(l) > 0){
    S <- S[,-l]
    F <- F[,-l] 
  }
  k <- rowSums(F)
  kn <- which(k == 1)
  if (length(kn) >0) {
    F <- F[-kn,]
  }
  b <- nrow(S)
  g <- colSums(S != 0)
  n <- colSums(S[,1:length(S)-1])
  l <- n/sum(n)
  fn <- g/b
  p <- which(l < 0.01 & fn[length(fn)-1]<=.5)
  if (length(p) >0) {
    F <- F[,-p]
    S <- S[,-p]
  }
  d <- colnames(S)
  d1<- which(d =="Chl.b")
  b <- rownames(F)
  b1<- which(b =="Chlorophytes")
  
  if(length(d1) ==0 & length(b1) >0){
    F <- F[-b1,]
  }
  c1 <- which(b =="Prasinophytes")
  
  if(length(d1) ==0 & length(c1) >0){
    F <- F[-c1,]
  }
  c <- rownames(F)
  c1 <- which(c =="Prasinophytes")
  d <- colnames(F)
  d1<- which(d =="Pra")
  if(length(d1) ==0 & length(c1) >0){
    F <- F[-c1,]
  }
  c <- rownames(F)
  c1 <- which(c =="Dinoflagellates-A")
  d <- colnames(F)
  d1<- which(d =="Per")
  if(length(d1) ==0 & length(c1) >0){
    F <- F[-c1,]
  }
  c <- rownames(F)
  c1 <- which(c =="Diatoms-A")
  d <- colnames(F)
  d1<- which(d =="Chl.c1")
  if(length(d1) ==0 & length(c1) >0){
    F <- F[-c1,]
  }
  c <- rownames(F)
  c1 <- which(c =="Dinoflagellates-A")
  d <- colnames(F)
  d1<- which(d =="Per")
  if(length(d1) ==0 & length(c1) >0){
    F <- F[-c1,]
  }
  c <- rownames(F)
  c1 <- which(c =="Syn")
  d <- colnames(F)
  d1<- which(d =="Zea")
  if(length(d1) ==0 & length(c1) >0){
    F <- F[-c1,]
  }
  
  c <- rownames(F)
  c1 <- which(c =="Cryptophytes")
  d <- colnames(F)
  d1<- which(d =="Allo")
  if(length(d1) ==0 & length(c1) >0){
    F <- F[-c1,]
  }
  
  c <- rownames(F)
  c1 <- which(c =="Haptophytes-H")
  d <- colnames(F)
  d1<- which(d =="X19hex")
  if(length(d1) ==0 & length(c1) >0){
    F <- F[-c1,]
  }
  c <- rownames(F)
  c1 <- which(c =="Haptophytes-L")
  d <- colnames(F)
  d1<- which(d =="X19hex")
  if(length(d1) ==0 & length(c1) >0){
    F <- F[-c1,]
  }  
  # c <- rownames(F)
  # c1 <- which(c =="Diatoms-A")
  # d <- colnames(F)
  # d1<- which(d =="Chl.c1")
  # if(length(d1) ==0 & length(c1) >0){
  #   F <- F[-c1,]
  # }
  
  
  
  d <- colnames(S)
  d1<- which(d =="X19but")
  b <- rownames(F)
  b1<- which(b =="Pelagophytes")
  if(length(d1) ==0 & length(b1) >0){
    F <- F[-b1,]
  }
  d <- colnames(S)
  d1<- which(d =="Chl.b")
  b <- rownames(F)
  b1<- which(b =="Prasinophytes")
  if(length(d1) == 0  & length(b1) >0){
    F <- F[-b1,]
  }
  k <- colSums(F)
  kn <- which(k == 0)
  if (length(kn) >0) {
    F <- F[,-kn]
    S <- S[,-kn]
  }
  
  
  return(list(S,F))
}
