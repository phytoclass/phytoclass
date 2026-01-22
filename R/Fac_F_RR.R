#' Part of the steepest descent algorithm that works to reduce error given 
#' the S and F matrices.
#' 
#' @keywords internal 
#'
#' @param Fmat A list containing `F matrix`, `RMSE` and `C matrix`
#' @param vary Indices of non-zero elements to vary in the optimization
#' @param S A matrix of samples (rows) and pigments (columns)
#' @param cm A vector of bounded weights for each pigment
#' @param fac_rr A numeric value (1, 2, or 3) to select which scaler values to use:
#'        1: (0.99, 1.01), 2: (0.98, 1.02), 3: (0.97, 1.03)
#' @param place A vector of all the indices of non-zero pigment ratios
#'
#' @return A list containing two elements:
#'   \code{1}: Updated F matrix after optimization
#'   \code{2}: Vector of indices
#'
#' @examples
#'  # Setup based on Minimise_elements_comb usage
#'  Fmat <- as.matrix(phytoclass::Fm)
#'  S <- as.matrix(phytoclass::Sm)
#'  cm <- as.numeric(phytoclass:::Bounded_weights(S))
#'  place <- which(Fmat[, seq(ncol(Fmat) - 2)] > 0)
#'
#'  # Get F.new from Conduit as done in Minimise_elements_comb
#'  f <- phytoclass:::Conduit(Fmat, place, S, cm, c_num = 3)
#'  F.new <- f[[1]]
#'
#'  # Set place1 as done in Minimise_elements_comb
#'  place1 <- place # place1 = place when c1_num != 1
#'
#'  # Run Fac_F_RR
#'  result <- phytoclass:::Fac_F_RR(F.new, vary = place, place = place1, S, cm, fac_rr = 3)
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
  
  # extract the pigments positions that, when adjusted, gives a better RMSE
  # ===
  # loops through each list of adjusted pigments and determines which adjustment
  # leads to a better RMSE based on the 4th list from `Replace_Rand` 
  cont  <- vapply(F.new, function(i) {i[[4]]}, logical(1))
  # extracts the index which has a better RMSE (i.e. TRUE)
  conts <- which(cont == 1) 
  
  # if the length is null, will either run through NNLS or recursively go to 
  # level one of Fac_F_RR
  # NOTE: this doesn't seem to be possible because when length(NULL) = 0 and is
  #       never null
  if (is.null(length(conts))) {
    if (fac_rr == 1) {
      F.new <- NNLS_MF(Fmat[[1]], S, cm)
      cont  <- vary
    } else {
      # TODO: remove because doesnt have enough inputs to work
      c1    <- Fac_F_RR(Fmat, place, fac_rr = 1)
      F.new <- c1[[1]]
      cont  <- c1[[2]]
    }
  }
  
  # when at least one pigment adjustment has a better RMSE, will replace all
  # pigment ratios in the Fmat with the better one and re-run NNLS
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
    # if no pigments return a better RMSE, depending on the number of iterations
    # in `Steepest_Descent`, will try a smaller range of scaler values
    if (fac_rr == 1) {
      # fac_rr = 1: smallest range, not re-running and extracting output only
      # from NNLS
      F.new <- NNLS_MF(Fmat[[1]], S, cm)
      cont  <- vary
    } else if (fac_rr == 2) {
      # fac_rr = 2: medium range going down to smallest range
      c1    <- Fac_F_RR(Fmat, place, S, cm, fac_rr = 1)
      F.new <- c1[[1]]
      cont  <- c1[[2]]
    } else if (fac_rr == 3) {
      # fac_rr = 3: largest range going down to smaller range
      c1    <- Fac_F_RR(Fmat, vary, place = place, S, cm, fac_rr = 2)
      F.new <- c1[[1]]
      cont  <- c1[[2]]  
    }
  }
  
  res <- list(F.new, cont)
  return(res)
}

# ============================================================================ #
# ---- old versions ---- #
# ============================================================================ #

#' #' Part of the steepest descent algorithm and work to reduce error given 
#' #' the S and F matrices.
#' #' 
#' #' @keywords internal 
#' #'
#' #' @param Fmat XX
#' #' @param vary XX
#' #' @param S XX
#' #' @param cm XX
#' #'
#' #' @return
#' #'
#' #' @examples
#' # Inputs are F and which elements to vary, should be all elements
#' Fac_F_RR1 <- function(Fmat, vary, S, cm){ 
#'   F.locs <- vector() 
#'   Fs <- Fmat[[1]]
#'   F.new <- lapply(vary, function(i){ # randomises every element in 'vary' 
#'     Replace_Rand(Fmat, i, S, cm, min.scaler = 0.99, max.scaler = 1.01)})
#'   # Shows which elements reduce error  
#'   cont <- lapply(1:length(F.new), function(i){ c <- which(length(F.new[[i]]) == 4)}) 
#'   conts <- which(cont==1)
#'   # Procedure for if no elements reduce error  
#'   if(!is.null(length(conts))){ 
#'     cont <- as.list(vary[conts])
#'     conts <- as.list(conts)
#'     Locs <-cont
#'     F.news <- sapply(conts, function(i) {
#'       sapply(cont, function(j){
#'         F.locs[[length(F.locs)+1]] <- F.new[[i]][[1]][[j]]
#'       })
#'     })
#'     if (length(F.news) > 0){
#'       if (length(F.news) >1 ){F.news <- diag(F.news)}
#'       else{F.news <- F.news[[1]]}      
#'       cont <- unlist(cont)
#'       F.new <- replace(Fmat[[1]],cont,F.news)
#'       F.new <- NNLS_MF(F.new, S, cm)
#'     }
#'     else{
#'       F.new <- NNLS_MF(Fmat[[1]], S, cm)
#'       cont <- vary
#'     }
#'   }
#'   else{
#'     F.new <- NNLS_MF(Fmat[[1]], S, cm)
#'     cont <- vary
#'   }
#'   res <- list(F.new,cont)
#'   return(res)
#'   
#' }
#' 
#' #' Part of the steepest descent algorithm and work to reduce error given 
#' #' the S and F matrices
#' #' 
#' #' @keywords internal
#' #'
#' #' @param Fmat   xx
#' #' @param vary    xx
#' #' @param place   xx
#' #' @param S   xx
#' #' @param cm   xx
#' #'
#' #' @return
#' #'
#' #' @examples
#' Fac_F_RR2 <- function(Fmat, vary, place, S, cm){
#'   F.locs <- vector()
#'   F.new <- lapply(vary, function(i){ 
#'     Replace_Rand(Fmat,i, S, cm, min.scaler = 0.98, max.scaler = 1.02)})
#'   cont <- lapply(1:length(F.new), function(i){ c <- which(length(F.new[[i]]) == 4)})
#'   conts <- which(cont==1)
#'   if(!is.null(length(conts))){
#'     cont <- as.list(vary[conts])
#'     conts <- as.list(conts)
#'     Locs <-cont
#'     F.news <- sapply(conts, function(i) {
#'       sapply(cont, function(j){
#'         F.locs[[length(F.locs)+1]] <- F.new[[i]][[1]][[j]]
#'       })
#'     })
#'     if (length(F.news) > 0){
#'       if (length(F.news) >1 ){F.news <- diag(F.news)}
#'       else{F.news <- F.news[[1]]} 
#'       cont <- unlist(cont)
#'       F.new <- replace(Fmat[[1]],cont,F.news)
#'       F.new <- NNLS_MF(F.new,S,cm)
#'     }
#'     else{
#'       C <- Fac_F_RR1(Fmat, place, S, cm)
#'       F.new <- C[[1]] 
#'       cont <- C[[2]]
#'     }
#'   }
#'   else{
#'     C <- Fac_F_RR1(Fmat, place, S, cm)
#'     F.new <- C[[1]] 
#'     cont <- C[[2]]
#'   }
#'   res <- list(F.new,cont)
#'   return(res)
#'   
#' }
#' 
#' #' Part of the steepest descent algorithm and work to reduce error given 
#' #' the S and F matrices
#' #' 
#' #' @keywords internal
#' #'
#' #' @param Fmat  xx
#' #' @param vary xx
#' #' @param place   xx   
#' #' @param S  xx
#' #' @param cm xx
#' #'
#' #' @return
#' #'
#' #' @examples
#' Fac_F_RR3 <- function(Fmat, vary, place, S, cm){
#'   F.locs <- vector()
#'   F.new <- lapply(vary, function(i) {
#'     Replace_Rand(Fmat, i, S, cm, min.scaler = 0.97, max.scaler = 1.03)
#'   })
#'   cont <- lapply(1:length(F.new), function(i){ c <- which(length(F.new[[i]]) == 4)})
#'   conts <- which(cont==1)
#'   if(!is.null(length(conts))){
#'     cont <- as.list(vary[conts])
#'     conts <- as.list(conts)
#'     Locs <-cont
#'     F.news <- sapply(conts, function(i) {
#'       sapply(cont, function(j){
#'         F.locs[[length(F.locs)+1]] <- F.new[[i]][[1]][[j]]
#'       })
#'     })
#'     if (length(F.news) > 0){
#'       if (length(F.news) >1 ){F.news <- diag(F.news)}
#'       else{F.news <- F.news[[1]]}
#'       cont <- unlist(cont)
#'       F.new <- replace(Fmat[[1]],cont,F.news)
#'       F.new <- NNLS_MF(F.new, S, cm)
#'     }
#'     else{
#'       C <- Fac_F_RR2(Fmat, vary, place, S, cm)
#'       F.new <- C[[1]] 
#'       cont <- C[[2]]
#'     }
#'   }
#'   else{
#'     C <- Fac_F_RR1(Fmat, place)
#'     F.new <- C[[1]] 
#'     cont <- C[[2]]
#'   }
#'   res <- list(F.new,cont)
#'   return(res)
#'   
#' }

