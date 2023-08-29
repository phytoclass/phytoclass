#' Turn each non-zero element of the F-matrix into a vector
#' 
#' @keywords internal
#'
#' @param F  xx
#'
#' @return
#'
#' @examples
vectorise <- function(F){
  g <- vector()
  for (i in F){
    if (i > -0){
      g[[length(g)+1]] <- i
      
    }
  }
  return(g)
}
