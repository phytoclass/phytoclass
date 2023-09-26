#' Performs the steepest descent algorithm for a set number of iterations
#' 
#' @keywords internal
#' 
#' @param Fmat xx
#' @param place xx
#' @param S   xx
#' @param cm   xx
#' @param num.loops   xx
#'
#' @return
#'
#' @examples
Steepest_Descent <- function(Fmat, place, S, cm, num.loops){ 
  loop <- 1
  F.new <- NNLS_MF(Fmat, S, cm)
  F.initial <- F.new
  for (i in 1:num.loops){ #should always be small. It would be nice to allow the 
    F.new <- Minimise_elements(F.initial[[1]], place, S, cm)
    loop = loop +1
    # print(loop)
    loop_2 <- 1
    # print(F.new[[2]])
    while (F.new[[2]] > F.initial[[2]]) {
      loop_2 = loop_2+1
      #print(loop_2)
      F.new <- Minimise_elements(F.initial[[1]], place, S, cm)
      if (loop_2 > 5){F.new <-  Minimise_elements1(F.initial[[1]], place,S,cm)} # If it doesn't work the first time, it randomises at a lower rate
      if (loop_2 > 10){F.new <- Minimise_elements2(F.initial[[1]], place,S,cm)} # and again, 
      
      #     if (loop_2 > 20){F.new <- Minimise_elements2(F.new[[1]])}
      if (loop_2 > 100){break} # it will continue for 100 itertions, and then stop
    }
    F.initial <- F.new
  }
  return(F.new)
}
