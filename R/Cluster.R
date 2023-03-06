#' Cluster things
#'
#' @param Data
#'
#' @return
#' @export
#'
#' @examples
Cluster <- function(Data) {
  standardise <- function(Data) {
    b <- Data
    b <- b[, 1:ncol(b) - 1]
    Chl <- Data[, ncol(Data)]
    b[b == 0] <- 1e-6
    b <- b / Chl
    v <- lapply(b, bestNormalize::boxcox)
    return(v)
  }

  v <- standardise(Data)

  L <- length(v)

  ndf <- list()
  for (i in 1:L) {
    ndf[[length(ndf) + 1]] <- data.frame(v[[i]][1])
  }

  S <- Data
  ndf <- do.call("cbind", ndf)
  colnames(ndf) <- colnames(S[, ncol(S) - 1])


  mscluster <- stats::dist(ndf, method = "manhattan")
  mv.hclust <- stats::hclust(mscluster, method = "ward.D2")

  plot(mv.hclust)

  ev.clust <- Data
  dynamicCut <- dynamicTreeCut::cutreeDynamic(mv.hclust,
    cutHeight = 70,
    minClusterSize = 12,
    method = "hybrid",
    distM = as.matrix(stats::dist(ndf, method = "manhattan")),
    deepSplit = 4,
    pamStage = TRUE,
    pamRespectsDendro = TRUE,
    useMedoids = FALSE,
    maxDistToLabel = NULL,
    maxPamDist = 50,
    respectSmallClusters = TRUE,
  )

  # NULL assignment to stop NOTE during the package "Check"
  #  -  no visible binding for global variable
  Clust <- NULL
  ev.clust$Clust <- dynamicCut

  L2 <- length(unique(ev.clust$Clust))
  L <- list()
  for (i in 1:L2) {
    L[[length(L) + 1]] <- dplyr::filter(ev.clust, Clust == i)
  }

  L

  e <- numeric()
  for (i in 1:length(L)) {
    e[[length(e) + 1]] <- length(L[[i]][[1]])
  }

  return(list(L, plot(mv.hclust)))
}
