#' Cluster things
#'
#' @param Data S (sample) matrix
#' @param min_cluster_size the minimum size required for a cluster
#'
#' @return A named list of length two. The first element "cluster.list"
#'  is a list of clusters, and the second element "cluster.plot" the 
#'  cluster analysis object (dendogram) that can be plotted.
#' @export
#'
#' @examples
#' Cluster.result <- Cluster(Sm, 14)
#' Cluster.result$cluster.list
#' plot(Cluster.result$cluster.plot)
Cluster <- function(Data, min_cluster_size) {
  
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


  mscluster <- stats::dist(ndf, method = "euclidean")
  mv.hclust <- stats::hclust(mscluster, method = "ward.D2")

  ev.clust <- Data
  # Change the minClusterSize argument below to change the number of clusters.
  # You might want to play around with it if you want ~ 6 clusters.
  # You could try setting it at 1/6th of the total sample number :)
  dynamicCut <- dynamicTreeCut::cutreeDynamic(mv.hclust,
    cutHeight = 70,
    minClusterSize = min_cluster_size,
    method = "hybrid",
    distM = as.matrix(stats::dist(ndf, method = "euclidean")), 
    deepSplit = 4,
    pamStage = TRUE, 
    pamRespectsDendro = TRUE,
    useMedoids = FALSE, 
    maxDistToLabel = NULL,
    maxPamDist = 50,
    respectSmallClusters = TRUE
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

  return(list(cluster.list = L, cluster.plot = mv.hclust))
}
