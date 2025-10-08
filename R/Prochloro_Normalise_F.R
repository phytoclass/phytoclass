#' Normalise F for prochloro
#'
#' @keywords internal
#'
#' @param Fmat 
#'
#' @return
#'
#' @examples
Prochloro_Normalise_F <- function(Fmat) {
  f_new <- as.matrix(Fmat)
  n <- nrow(f_new)
  p <- ncol(f_new)
  
  # Identify Pro row (prefer rowname; else assume last row)
  i_pro <- which(tolower(rownames(f_new)) %in% c("pro", "prochlorococcus", "prochlorococcus-1", "pro-1"))
  if (length(i_pro) != 1) i_pro <- n
  i_nonpro <- setdiff(seq_len(n), i_pro)
  
  chla   <- f_new[, p]        # last column is Chl a (Tchla)
  dvchla <- f_new[, p - 1]    # second last is dvChl a
  
  # --- Row-wise scaling so biomass pigment == 1 ---
  # Non-Pro groups: scale by Chl a of that row
  if (length(i_nonpro) > 0) {
    denom <- chla[i_nonpro]
    denom[!is.finite(denom) | denom == 0] <- 1   # guard
    f_new[i_nonpro, ] <- f_new[i_nonpro, , drop = FALSE] / denom
  }
  
  # Pro row: scale by its dvChl a
  dv_pro <- dvchla[i_pro]
  if (!is.finite(dv_pro) || dv_pro == 0) dv_pro <- 1
  f_new[i_pro, ] <- f_new[i_pro, , drop = FALSE] / dv_pro
  
  # Enforce exact biomass markers after scaling
  f_new[i_pro, p - 1] <- 1      # Pro: dvChl a == 1
  # (Chl a for Pro should already be 0; keep it as-is.)
  
  # Final step mirrors your original pipeline:
  # compute row sums AFTER pigment scaling; return the *pre-row-sum* scaled version via Fn <- Fn * F.sum
  f_sum <- rowSums(f_new)
  f_norm <- f_new / f_sum
  return(list(as.matrix(f_norm), f_sum))
}
