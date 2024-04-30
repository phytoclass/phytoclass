#' Remove any column values that average 0. Further to this, also remove
#' phytoplankton groups from the F matrix if their diagnostic pigment
#' isn’t present. 
#'
#' @param S   Sample data matrix – a matrix of pigment samples
#' @param Fmat   Pigment to Chl a matrix
#'
#' @return Named list with new S and Fmat matrices
#' @export
#'
#' @examples
#' MC <- Matrix_checks(Sm, Fm)  
#' Snew <- MC$Snew
# 'Fnew <- MC$Fnew 
#' 

Matrix_checks <- function(S, Fmat){
  # Only keep columns of Fmat that are in S
  S.colnames <- colnames(S)
  Fmat.colnames <- colnames(Fmat)
  keep.these.columns <- (Fmat.colnames %in% S.colnames)
  Fmat[, keep.these.columns]
  #
  ba <- rownames(Fmat)
  ba1<- which(ba =="Syn")
  if (ncol(S) > ncol(Fmat)){
    S <- subset(S, select = c(colnames(Fmat)))}
  b <- colSums(S)
  g <- mean(S[,ncol(S)])
  ba <- rownames(Fmat)
  ba1<- which(ba =="Syn")
  
  b <- nrow(S)
  c <- which(b ==0)
  g <- colSums(S != 0)
  l <- which(g/b <=.01)
  if(length(l) > 0){
    S <- S[,-l]
    Fmat <- Fmat[,-l] 
  }
  k <- rowSums(Fmat)
  kn <- which(k == 1)
  if (length(kn) >0) {
    Fmat <- Fmat[-kn,]
  }
  b <- nrow(S)
  g <- colSums(S != 0)
  n <- colSums(S[,1:ncol(S)-1])
  l <- n/sum(n)
  fn <- g/b
  p <- which(l < 0.01  & fn[1:length(fn) - 1] <= 0.5)
  #if (length(p) >0) {
  #  F <- F[,-p]
  #  S <- S[,-p]
  #}
  d <- colnames(S)
  d1<- which(d =="Chl_b")
  b <- rownames(Fmat)
  b1<- which(b =="Chlorophytes")
  
  if(length(d1) ==0 & length(b1) >0){
    Fmat <- Fmat[-b1,]
  }
  c1 <- which(b =="Prasinophytes")
  
  if(length(d1) ==0 & length(c1) >0){
    Fmat <- Fmat[-c1,]
  }
  c <- rownames(Fmat)
  c1 <- which(c =="Prasinophytes")
  d <- colnames(Fmat)
  d1<- which(d =="Pra")
  if(length(d1) ==0 & length(c1) >0){
    Fmat <- Fmat[-c1,]
  }
  c <- rownames(Fmat)
  c1 <- which(c =="Dinoflagellates-1")
  d <- colnames(Fmat)
  d1<- which(d =="Per")
  if(length(d1) ==0 & length(c1) >0){
    Fmat <- Fmat[-c1,]
  }
  c <- rownames(Fmat)
  c1 <- which(c =="Diatoms-1")
  d <- colnames(Fmat)
  d1<- which(d =="Chl.c1")
  if(length(d1) ==0 & length(c1) >0){
    Fmat <- Fmat[-c1,]
  }
  c <- rownames(Fmat)
  c1 <- which(c =="Diatoms-2")
  d <- colnames(Fmat)
  d1<- which(d =="Fuco")
  if(length(d1) ==0 & length(c1) >0){
    Fmat <- Fmat[-c1,]
  }
  c <- rownames(Fmat)
  c1 <- which(c =="Dinoflagellates-1")
  d <- colnames(Fmat)
  d1<- which(d =="Per")
  if(length(d1) ==0 & length(c1) >0){
    Fmat <- Fmat[-c1,]
  }
  c <- rownames(Fmat)
  c1 <- which(c =="Syn")
  d <- colnames(Fmat)
  d1<- which(d =="Zea")
  if(length(d1) ==0 & length(c1) >0){
    Fmat <- Fmat[-c1,]
  }
  
  c <- rownames(Fmat)
  c1 <- which(c =="Cryptophytes")
  d <- colnames(Fmat)
  d1<- which(d =="Allo")
  if(length(d1) ==0 & length(c1) >0){
    Fmat <- Fmat[-c1,]
  }
  
  c <- rownames(Fmat)
  c1 <- which(c =="Haptophytes-H")
  d <- colnames(Fmat)
  d1<- which(d =="X19hex")
  if(length(d1) ==0 & length(c1) >0){
    Fmat <- Fmat[-c1,]
  }
  c <- rownames(Fmat)
  c1 <- which(c =="Haptophytes-L")
  d <- colnames(Fmat)
  d1<- which(d =="X19hex")
  if(length(d1) ==0 & length(c1) >0){
    Fmat <- Fmat[-c1,]
  }  
  c <- rownames(Fmat)
  c1 <- which(c =="Diatoms-1")
  d <- colnames(Fmat)
  d1<- which(d =="Fuco")
  if(length(d1) ==0 & length(c1) >0){
    Fmat <- Fmat[-c1,]
  }
  
  d <- colnames(S)
  d1<- which(d =="X19but")
  b <- rownames(Fmat)
  b1<- which(b =="Pelagophytes")
  if(length(d1) ==0 & length(b1) >0){
    Fmat <- Fmat[-b1,]
  }
  d <- colnames(S)
  d1<- which(d =="Chl.b")
  b <- rownames(Fmat)
  b1<- which(b =="Prasinophytes")
  if(length(d1) == 0  & length(b1) >0){
    Fmat <- Fmat[-b1,]
  }
  k <- colSums(Fmat)
  kn <- which(k == 0)
  if (length(kn) >0) {
    Fmat <- Fmat[,-kn]
    S <- S[,-kn]
  }
  return(list(Snew = as.matrix(S), Fnew = as.matrix(Fmat)))
}
