#' Part of the steepest descent algorithm and work to reduce error given 
#' the S and F matrices.
#' 
#' @keywords internal 
#'
#' @param Fmat A list containing `F matrix`, `RMSE` and `C matrix`
#' @param vary Are the index of the non-zero integers to vary
#' @param S A matrix of samples (rows) and pigments (columns)
#' @param cm A vector of bounded weights for each pigment
#' @param fac_rr A numeric to call which scaler values to use
#' @param place A vector of all the indexs of non-zero pigment ratios
#'
#' @return
#'
#' @examples
# Inputs are F and which elements to vary, should be all elements
Fac_F_RR <- function(Fmat, vary, S, cm, fac_rr = c(1, 2, 3), place = NULL) {
  
  # which min and max scale values to use
  rand_var <-
    switch(
      fac_rr,
      list(min.scaler = 0.99, max.scaler = 1.01),
      list(min.scaler = 0.98, max.scaler = 1.02),
      list(min.scaler = 0.97, max.scaler = 1.03)
    )
  
  # randomises every element in 'vary' then retest RMSE to prev RMSE
  F.new <- lapply(
    vary, 
    function(i) {
      Replace_Rand(
        Fmat, i, S, cm,
        min.scaler = rand_var$min.scaler,
        max.scaler = rand_var$max.scaler
      )
    }
  )
  
  cont  <- vapply(F.new, function(i) {i[[4]]}, logical(1))
  conts <- which(cont == 1)
  
  if (is.null(length(conts))) {
    if (fac_rr == 1) {
      F.new <- NNLS_MF(Fmat[[1]], S, cm)
      cont  <- vary
    } else {
      # TODO: remove because doesnt have enough inputs to work
      C1    <- Fac_F_RR(Fmat, place, fac_rr = 1)
      F.new <- C1[[1]]
      cont  <- C1[[2]]
    }
  }
  
  if (length(conts) > 0) {
    F.news <- vector(length = length(conts)) # initialize new ratios
    for (i in 1:length(conts)) {
      pig_ind     <- vary[conts][i]
      F_ind       <- conts[i]
      F.news[[i]] <- F.new[[F_ind]][[1]][pig_ind]
    }
    F.new <- replace(Fmat[[1]], vary[conts], F.news)
    F.new <- NNLS_MF(F.new, S, cm)
  } else {
    
    if (fac_rr == 1) {
      F.new <- NNLS_MF(Fmat[[1]], S, cm)
      cont  <- vary
    } else if (fac_rr == 2) {
      
      C1    <- Fac_F_RR(Fmat, place, S, cm, fac_rr = 1)
      F.new <- C1[[1]]
      cont  <- C[[2]]
    } else if (fac_rr == 3) {
      
      C1    <- Fac_F_RR(Fmat, vary, place = place, S, cm, fac_rr = 2)
      F.new <- C1[[1]]
      cont  <- C1[[2]]  
    }
  }
  
  res <- list(F.new, cont)
  return(res)
}