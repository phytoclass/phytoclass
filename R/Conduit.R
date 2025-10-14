#' Conduit between minimise_elements function and Fac_F_R 
#' of steepest descent algorithm.
#' 
#' @keywords internal
#'
#' @param Fmat  xx   
#' @param place  xx
#' @param S  xx
#' @param cm  xx
#' @return
#'
#' @examples
Conduit <- function(Fmat, place, S, cm, c_num = c(1, 2, 3)) {
  # run NNLS to get previous RMSE
  F.old <- NNLS_MF(Fmat, S, cm)
  
  # 
  place1 <- NULL
  if (c_num != 1) {place1 <- place }
  
  F.news <- Fac_F_RR(F.old, vary = place, place = place1, S, cm, fac_rr = c_num)
  
  res <- list(F.news[[1]], F.news[[2]], F.old)
  return(res)
}

# ============================================================================ #
# ---- old versions ---- #
# ============================================================================ #

#' #' Conduit between minimise_elements function and Fac_F_R 
#' #' of steepest descent algorithm.
#' #' 
#' #' @keywords internal
#' #'
#' #' @param Fmat  xx   
#' #' @param place  xx
#' #' @param S  xx
#' #' @param cm  xx
#' #' @return
#' #'
#' #' @examples
#' Conduit_1 <- function(Fmat, place, S, cm){
#'   F.locs <- vector()
#'   F.old <- NNLS_MF(Fmat, S, cm)
#'   F.news <- Fac_F_RR1(F.old, place, S, cm)
#'   F.new <- F.news[[1]]
#'   n <- F.news[[2]]
#'   res <- list(F.new,n, F.old)
#'   return(res)
#' }
#' 
#' #' Conduit between minimise_elements function and Fac_F_R 
#' #' of steepest descent algorithm.
#' #' 
#' #' @keywords internal
#' #' 
#' #' @param Fmat    xx 
#' #' @param place   xx
#' #' @param S   xx
#' #' @param cm   xx
#' #'
#' #' @return
#' #'
#' #' @examples
#' Conduit_2 <- function(Fmat, place, S, cm){
#'   F.locs <- vector()
#'   F.old <- NNLS_MF(Fmat, S, cm)
#'   F.news <- Fac_F_RR2(F.old, vary = place, place, S, cm)
#'   F.new <- F.news[[1]]
#'   n <- F.news[[2]]
#'   res <- list(F.new, n, F.old)
#'   return(res)
#' }
#' 
#' #' Conduit between minimise_elements function and Fac_F_R 
#' #' of steepest descent algorithm.
#' #' 
#' #' @keywords internal
#' #' 
#' #' @param Fmat   xx
#' #' @param place xx
#' #' @param S  xx
#' #' @param cm xx
#' #'
#' #' @return
#' #'
#' #' @examples
#' Conduit_3 <- function(Fmat, place, S, cm){
#'   F.old <- NNLS_MF(Fmat, S, cm)
#'   F.news <- Fac_F_RR3(F.old, vary = place, place, S, cm)
#'   F.new <- F.news[[1]]
#'   n <- F.news[[2]]
#'   res <- list(F.new, n, F.old)
#'   return(res)
#' }
