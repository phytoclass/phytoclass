#' Sets the default minimum and maximum values for phytoplankton groups 
#' pigment ratios. To use this function, pigment and phytoplankton group
#' names will need to fit the naming criteria of phytoclass. 
#'  
#' @keywords internal
#' 
#' @param min_max A data.frame with 3 4 columns for class, pigment, min and max.
#' @param Fmat    F matrix
#'
#' @return
#'
#' @examples
#' 
Default_min_max <- function(min_max, Fmat) {
  
  # selects all non-zero pigment ratios and indexes the phyto names and pigments
  k     <- which(Fmat > 0, TRUE)
  RName <- rownames(Fmat)[k[, 1]]
  CName <- colnames(Fmat)[k[, 2]]

  # Find indices in min_max matching each taxa-pigment pair;
  # throw error if any pair is missing in min_max
  vecs    <- vector(length = nrow(k))

  for (i in seq_len(nrow(k))) {
    idx <- which(RName[i] == min_max[, 1] & CName[i] == min_max[, 2])
    
    # add -1 if missing phyto - pigment pair
    if (length(idx) == 0) {
      idx <- -1
    }
    vecs[i] <- idx
  }

  # collect all missing pairs and error
  if (any(vecs == -1)) {
    mis_idx <- which(vecs == -1)
    stop(
      call. = FALSE,
      c(
        "`Default_min_max()` - Your F matrix includes a missing taxa-pigment ", 
        "pair in the `min_max` matrix.\n",
        "Missing pairs:\n",
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
