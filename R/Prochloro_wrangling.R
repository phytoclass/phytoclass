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
# Prochloro_Wrangling <- function (Fl, min.val, max.val) 
# {
#   Fd             <- Fl
#   Fmin           <- as.matrix(Fd)
#   Fmin           <- Fmin[, 1:ncol(Fmin) - 1]
#   Fmin[Fmin > 0] <- min.val
#   Fmax           <- as.matrix(Fd)
#   Fmax           <- Fmax[, 1:ncol(Fmax) - 1]
#   Fmax[Fmax > 0] <- max.val
#   chlv           <- Fd[, ncol(Fd)]
#   chlvp          <- Fd[, ncol(Fd)-1]
#   chlep          <- Fd[, ncol(Fd)-1]
#   chlep[length(chlep)] <- 1
#   
#   Fmin[1:nrow(Fmin)-1,1:ncol(Fmin)-1] <- Fmin[1:nrow(Fmin)-1,1:ncol(Fmin)-1] * chlv[1:length(chlv)-1]
#   Fmin[nrow(Fmin),1:ncol(Fmin)] <- Fmin[nrow(Fmin),1:ncol(Fmin)] * chlvp[length(chlvp)]
#   Fmin <- cbind(Fmin[,1:ncol(Fmin)-1], chlep, chlv)
#   
#   Fmax[1:nrow(Fmax)-1,1:ncol(Fmax)-1] <- Fmax[1:nrow(Fmax)-1,1:ncol(Fmax)-1] * chlv[1:length(chlv)-1]
#   Fmax[nrow(Fmax),1:ncol(Fmax)] <- Fmax[nrow(Fmax),1:ncol(Fmax)] * chlvp[length(chlvp)]
#   Fmax <- cbind(Fmax[,1:ncol(Fmax)-1],chlep, chlv)
#  
#   Fmin <- vectorise(Fmin[, 1:ncol(Fmin) - 1])
#   Fmax <- vectorise(Fmax[, 1:ncol(Fmax) - 1])
#   SE   <- vectorise(Fd[, 1:ncol(Fd) - 1])
#   res  <- list(Fmin, Fmax, SE, chlv,chlvp)
# }

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

