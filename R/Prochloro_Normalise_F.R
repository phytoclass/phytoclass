#' Normalise F for prochloro
#'
#' @keywords internal
#'
#' @param Fmat 
#'
#' @return
#'
#' @examples
Prochloro_Normalise_F <- function (Fmat) {
  F <- as.matrix(Fmat)
  n <- nrow(F); p <- ncol(F)
  
  # Identify Pro row (prefer rowname; else assume last row)
  i_pro <- which(tolower(rownames(F)) %in% c("pro", "prochlorococcus", "prochlorococcus-1", "pro-1"))
  if (length(i_pro) != 1) i_pro <- n
  i_nonpro <- setdiff(seq_len(n), i_pro)
  
  chla   <- F[, p]        # last column is Chl a (Tchla)
  dvchla <- F[, p - 1]    # second last is dvChl a
  
  # --- Row-wise scaling so biomass pigment == 1 ---
  # Non-Pro groups: scale by Chl a of that row
  if (length(i_nonpro) > 0) {
    denom <- chla[i_nonpro]
    denom[!is.finite(denom) | denom == 0] <- 1   # guard
    F[i_nonpro, ] <- F[i_nonpro, , drop = FALSE] / denom
  }
  
  # Pro row: scale by its dvChl a
  dv_pro <- dvchla[i_pro]
  if (!is.finite(dv_pro) || dv_pro == 0) dv_pro <- 1
  F[i_pro, ] <- F[i_pro, , drop = FALSE] / dv_pro
  
  # Enforce exact biomass markers after scaling
  F[i_pro, p - 1] <- 1      # Pro: dvChl a == 1
  # (Chl a for Pro should already be 0; keep it as-is.)
  
  # Final step mirrors your original pipeline:
  # compute row sums AFTER pigment scaling; return the *pre-row-sum* scaled version via Fn <- Fn * F.sum
  F.sum <- rowSums(F)
  F.norm <- F / F.sum
  list(as.matrix(F.norm), F.sum)
}
