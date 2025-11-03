#' Sets the default minimum and maximum values for phytoplankton groups 
#' pigment ratios. To use this function, pigment and phytoplankton group
#' names will need to fit the naming criteria of phytoclass. 
#'  
#' @keywords internal
#' 
#' @param min_max A data.frame with 4 columns for class, pigment, min and max values
#' @param Fmat F matrix with phytoplankton groups as rows and pigments as columns
#'
#' @return A list containing two elements:
#'   \item{[[1]]}{Vector of minimum values for each non-zero pigment ratio}
#'   \item{[[2]]}{Vector of maximum values for each non-zero pigment ratio}
#'
#' @examples
#' # Create a sample F matrix
#' Fmat <- matrix(c(0.5, 0, 0.3,
#'                  0, 0.4, 0.2), 
#'                nrow=2, byrow=TRUE,
#'                dimnames=list(c("diatoms", "cyano"),
#'                             c("chl_a", "chl_b", "chl_c")))
#' 
#' # Create min_max data frame
#' min_max <- data.frame(
#'   class=c("diatoms", "diatoms", "cyano", "cyano"),
#'   pigment=c("chl_a", "chl_c", "chl_b", "chl_c"),
#'   min=c(0.4, 0.2, 0.3, 0.1),
#'   max=c(0.6, 0.4, 0.5, 0.3)
#' )
#' 
#' result <- Default_min_max(min_max, Fmat)
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
