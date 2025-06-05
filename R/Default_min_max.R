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
#' 
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

Default_min_max <- function(min_max, Fmat) {
  
  # selects all non-zero pigment ratios and indexes the phyto names and pigments
  k <- which(Fmat > 0, TRUE)
  RName <- rownames(Fmat)[k[, 1]]
  CName <- colnames(Fmat)[k[, 2]]

  # Find indices in min_max matching each taxa-pigment pair;
  # throw error if any pair is missing in min_max
  vecs    <- vector(length = nrow(k))

  for (i in seq_len(nrow(k))) {
    idx <- which(RName[i] == min_max[, 1] & CName[i] == min_max[, 2])
    if (length(idx) == 0) {
      idx <- -1
    }
    vecs[i] <- idx
  }

  if (any(vecs == -1)) {
    mis_idx <- which(vecs == -1)
    stop(
      call. = FALSE,
      c(
        "Your F matrix includes a missing taxa-pigment pair in the `min_max` matrix.\n",
        sprintf(
          "%-20.20s : %7.7s\n",
          RName[mis_idx],
          CName[mis_idx]
        )
      )
    )
  }

  min <- min_max[vecs, 3]
  max <- min_max[vecs, 4]

  return(list(min, max))
}
