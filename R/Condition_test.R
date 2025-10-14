#' Calculate the condition number ...
#' 
#' @keywords internal
#'
#' @param S XX     
#' @param Fn  XX   
#' @param min.val XX     
#' @param max.val XX       
#'
#' @return
#'
#' @examples
Condition_test <- function(S, Fn, min.val = NULL, max.val = NULL) {
  if (is.null(min.val) & is.null(max.val)) {
    min_max <- Default_min_max(phytoclass::min_max, Fn)
    min.val <- min_max[[1]]
    max.val <- min_max[[2]]
  }
  condition_number <- function(S, f_mat, min.val, max.val) {
    f_non_zero <- length(vectorise(f_mat))
    rand       <- vector(length = f_non_zero)
    for (i in seq(f_non_zero)) {
      rand[i] <- stats::runif(1, min = min.val[i], max = max.val[i])
    }
    f_mat[f_mat > 0] <- rand
    return(kappa(f_mat %*% t(S)))
  }
  sn <- replicate(n = 1000, condition_number(S, Fn, min.val, max.val))
  return(mean(sn))
}
