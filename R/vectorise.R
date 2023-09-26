#' Turn each non-zero element of the F-matrix into a vector
#' 
#' @keywords internal
#'
#' @param Fmat  xx
#'
#' @return
#'
#' @examples
vectorise <- function(Fmat){
  g <- vector()
  for (i in Fmat){
    if (i > -0){
      g[[length(g)+1]] <- i
      
    }
  }
  return(g)
}
