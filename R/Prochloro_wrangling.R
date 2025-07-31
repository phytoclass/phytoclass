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
Prochloro_Wrangling <- function (Fl, min.val, max.val) 
{
  Fd <- Fl
  Fmin <- as.matrix(Fd)
  Fmin <- Fmin[, 1:ncol(Fmin) - 1]
  Fmin[Fmin > 0] <- min.val
  Fmax <- as.matrix(Fd)
  Fmax <- Fmax[, 1:ncol(Fmax) - 1]
  Fmax[Fmax > 0] <- max.val
  chlv <- Fd[, ncol(Fd)]
  chlvp <- Fd[, ncol(Fd)-1]
  chlep <- Fd[, ncol(Fd)-1]
  chlep[length(chlep)] <- 1
  Fmin[1:nrow(Fmin)-1,1:ncol(Fmin)-1] <- Fmin[1:nrow(Fmin)-1,1:ncol(Fmin)-1] * chlv[1:length(chlv)-1]
  Fmin[nrow(Fmin),1:ncol(Fmin)] <- Fmin[nrow(Fmin),1:ncol(Fmin)] * chlvp[length(chlvp)]
  Fmin <- cbind(Fmin[,1:ncol(Fmin)-1], chlep, chlv)
  Fmax[1:nrow(Fmax)-1,1:ncol(Fmax)-1] <- Fmax[1:nrow(Fmax)-1,1:ncol(Fmax)-1] * chlv[1:length(chlv)-1]
  Fmax[nrow(Fmax),1:ncol(Fmax)] <- Fmax[nrow(Fmax),1:ncol(Fmax)] * chlvp[length(chlvp)]
  Fmax <- cbind(Fmax[,1:ncol(Fmax)-1],chlep, chlv)
  Fmin <- vectorise(Fmin[, 1:ncol(Fmin) - 1])
  Fmax <- vectorise(Fmax[, 1:ncol(Fmax) - 1])
  SE <- vectorise(Fd[, 1:ncol(Fd) - 1])
  res <- list(Fmin, Fmax, SE, chlv,chlvp)
}


