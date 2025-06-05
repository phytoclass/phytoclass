#' Sets the default minimum and maximum values for phytoplankton groups 
#' pigment ratios. To use this function, pigment and phytoplankton group
#' names will need to fit the naming criteria of phytoclass. 
#'  
#' @keywords internal
#' 
#' @param min_max   xx
#' @param Fmat    xx
#' @param place       xx
#'
#' @return
#'
#' @examples


Default_min_max <- function(min_max, Fmat, place){
  k <- list()
  for (i in place){
    k[[length(k)+1]] <- arrayInd(i, dim(Fmat))
  }
  
  RName <- vector()
  CName <- vector()
  for (i in 1:length(k)){
    RName[[length(RName)+1]] <- rownames(Fmat)[k[[i]][,1]]
    CName[[length(CName)+1]] <- colnames(Fmat)[k[[i]][,2]]
  }
  
  
  # Find indices in min_max matching each taxa-pigment pair;
  # throw error if any pair is missing in min_max
  vecs <- vector()
  for (i in 1:length(k)){
    idx <- which(RName[[i]] == min_max[,1] & CName[[i]] == min_max[,2])
    if (length(idx) == 0) {
      stop(paste0("Your F matrix includes an unexpected taxa-pigment pair for ", 
                  RName[[i]], " - ", CName[[i]], 
                  ". This pair is not in the min_max matrix."))
    }
    vecs[[length(vecs)+1]] <- idx
  }
  
  min <- min_max[vecs,3]
  max <- min_max[vecs,4]
  
  return(list(min,max))
}
