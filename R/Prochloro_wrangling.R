#' Prochloro Wrangling
#' 
#' @keywords internal
#'
#' @param Fl 
#' @param min.val 
#' @param max.val 
#'
#' @return
#'
#' @examples
Prochloro_Wrangling <- function(Fl, min.val, max.val) {

  # set up initial F, Fmin, Fmax, matrix by removing Tchla column
  Fd <- Fl
  Fmin <- Fmax <- as.matrix(Fl)[, -ncol(Fl)]

  Fmin[Fmin > 0] <- min.val
  Fmax[Fmax > 0] <- max.val

  chlv  <- Fl[, ncol(Fl)]      # Chl a
  chlvp <- Fl[, ncol(Fl) - 1]  # dvChl a
  chlep <- Fl[, ncol(Fl) - 1]
  chlep[length(chlep)] <- 1

  # Non-Prochloro rows scaled by Chl a (exclude dvChl column)
  Fmin[-nrow(Fmin), -ncol(Fmin)] <- Fmin[-nrow(Fmin), -ncol(Fmin)] * chlv[-length(chlv)]
  Fmax[-nrow(Fmax), -ncol(Fmax)] <- Fmax[-nrow(Fmax), -ncol(Fmax)] * chlv[-length(chlv)]

  # Prochlorococcus row scaled by dvChl a
  Fmin[nrow(Fmin), ] <- Fmin[nrow(Fmin), ] * chlvp[length(chlvp)]
  Fmax[nrow(Fmax), ] <- Fmax[nrow(Fmax), ] * chlvp[length(chlvp)]

  # Vectorised CORE (base-R column-major) and remove zeros
  Fmin_vec <- vectorise(Fmin)
  Fmax_vec <- vectorise(Fmax)
  SE_vec   <- vectorise(Fl[, 1:(ncol(Fmin) - 2)])
  
  return(list(Fmin_vec, Fmax_vec, SE_vec, chlv, chlvp))
}
