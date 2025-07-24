#' Remove any column values that average 0. Further to this, also remove
#' phytoplankton groups from the F matrix if their diagnostic pigment
#' isn’t present. 
#'
#' @param S   Sample data matrix – a matrix of pigment samples
#' @param Fmat   Pigment to Chl a matrix
#'
#' @return Named list with new S and Fmat matrices
#' @export
#'
#' @examples
#' MC <- Matrix_checks(Sm, Fm)  
#' Snew <- MC$Snew
# 'Fnew <- MC$Fnew 
#' 
Matrix_checks <- function(S, Fmat) {
  
  # check if non-numeric and if all groups have more than 1 pigment
  f_rowsum <- tryCatch({
    rowSums(Fmat)
  }, error = function(e) {
    if (grepl("'x' must be numeric", e$message)) {
      stop("Error in rowSums(Fmat) : 'x' must be numeric. \n\n The F-Matrix contains non-numeric data.")
    } else {
      stop(e)
    }
  })
  
  f_mat <- Fmat[!(f_rowsum <= 1),] # remove groups with only 1 pigment
  f_mat <- f_mat[, colSums(Fmat) > 0] # remove columns if all are 0
  
  # if column names dont match, select the columns that overlap with Fmat
  # also reorders to match Fmat
  s_mat <- S[,intersect(colnames(f_mat), colnames(S))]
  f_mat <- f_mat[,intersect(colnames(f_mat), colnames(S))]
  
  # remove pigments when occur in less than 1% of samples 
  s_nrow   <- nrow(s_mat)
  s_colsum <- colSums(s_mat != 0) # number of present pigments
  s_indx   <- which(s_colsum / s_nrow <= -0.01) # less than 1% have value for this pigments?
  s_indx   <- !(s_colsum / s_nrow <= 0.01) # less than 1% have value for this pigments?
  keep_col <- colnames(s_mat)[!(s_colsum / s_nrow <= 0.01)] # pigments to keep
  s_mat    <- s_mat[,keep_col]
  f_mat    <- f_mat[,keep_col]
  
  
  # check S matrix that pigments occur in > 50% of samples and total concentration
  # is >= 1% of total pigment pool
  # NOTE: removal of pigments is not implemented
  pig_percent_gt0 <- colSums(s_mat != 0)[-ncol(s_mat)] / s_nrow # pigment amount percent greater than 0 
  pig_conc        <- colSums(s_mat[,-ncol(s_mat)]) # sums all non-chlor-a pigments
  pig_percent     <- pig_conc / sum(pig_conc) # percent of total pigment conc
  # f_mat <- f_mat[, !(pig_percent < 0.01  & pig_percent_gt0 <= 0.5)]
  # s_mat <- s_mat[, !(pig_percent < 0.01  & pig_percent_gt0 <= 0.5)]
  
  check_mat <-
    matrix(
      c(
        c(
          "Chlorophytes", "Prasinophytes", "Prasinophytes", "Dinoflagellates-1",
          "Diatoms-1", "Diatoms-2", "Syn", "Cryptophytes", "Haptophytes-H",
          "Haptophytes-L", "Diatoms-1", "Pelagophytes", "Prasinophytes"
        ),
        c(
          "Chl_b", "Chl_b", "Pra", "Per", "Chl_c1", "Fuco", "Zea", "Allo", 
          "X19hex", "X19hex", "Fuco", "X19but", "Chl_b"
        )
      ),
      ncol = 2
    )
  
  for (i in seq(nrow(check_mat))) {
    phyto_row <- which(rownames(f_mat) == check_mat[i, 1])
    pig_col   <- which(colnames(s_mat) == check_mat[i, 2])
    if (length(pig_col) == 0 & length(phyto_row) > 0) {
      print("removing something")
      f_mat <- f_mat[-phyto_row,]
    }
  }
  
  # final check to remove F cols with no pigments
  kn <- which(colSums(f_mat) == 0)
  if (length(kn) > 0) {
    f_mat <- f_mat[,-kn]
    s_mat <- s_mat[,-kn]
  }
  
  return(list(Snew = as.matrix(s_mat), Fnew = as.matrix(f_mat)))
  
}

