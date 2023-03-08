#' Sets the default minimum and maximum values for phytoplankton groups 
#' pigment ratios. To use this function, pigment and phytoplankton group
#' names will need to fit the naming criteria of phytoclass. 
#'  
#' 
#' @param min_max   
#' @param F    
#' @param place      
#'
#' @return
#' @export
#'
#' @examples


Default_min_max <- function(min_max, F, place){
  k <- list()
  for (i in place){
    k[[length(k)+1]] <- arrayInd(i, dim(F))
  }
  
  RName <- vector()
  CName <- vector()
  for (i in 1:length(k)){
    RName[[length(RName)+1]] <- rownames(F)[k[[i]][,1]]
    CName[[length(CName)+1]] <- colnames(F)[k[[i]][,2]]
  }
  
  vecs <- vector()
  for (i in 1:length(k)){
    vecs[[length(vecs)+1]]<-which(RName[[i]] == min_max[,1] & CName[[i]] == min_max[,2])
  }
  
  min <- min_max[vecs,3]
  max <- min_max[vecs,4]
  
  return(list(min,max))
}
