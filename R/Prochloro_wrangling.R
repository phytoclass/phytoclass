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

  Fmin <- Fmax   <- Fl[, -ncol(Fl)]
  Fmin[Fmin > 0] <- min.val
  Fmax[Fmax > 0] <- max.val
  chlv           <- Fl[, ncol(Fl)] # chla
  chlvp <- chlep <- Fl[, ncol(Fl) - 1] # 2nd to last col dvchla
  chlep[length(chlep)] <- 1

  F_row <- nrow(Fmin)
  F_col <- ncol(Fmin)
  
  Fmin[-F_row, -F_col] <- Fmin[-F_row, -F_col] * chlv[-length(chlv)]
  Fmax[-F_row, -F_col] <- Fmax[-F_row, -F_col] * chlv[-length(chlv)]
  
  Fmin[F_row, ] <- Fmin[F_row, ] * chlvp[length(chlvp)]
  Fmax[F_row, ] <- Fmax[F_row, ] * chlvp[length(chlvp)]
  
  Fmin <- vectorise(cbind(Fmin[, -F_col], chlep))
  Fmax <- vectorise(cbind(Fmax[, -F_col], chlep))
  SE   <- vectorise(Fl[, -ncol(Fl)])

  return(list(Fmin, Fmax, SE, chlv,chlvp))
}

