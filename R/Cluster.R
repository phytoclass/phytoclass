#' Cluster things
#'
#' @param Data S (sample) matrix
#' @param min_cluster_size the minimum number of samples required for a cluster
#' @param row_ids A vector of custom row names to be added to dendrogram
#' @param dist_method Distance metric to be used in `stats::dist`. This should 
#'                    be one of "euclidean", "maximum", "manhattan", "canberra", 
#'                    "binary" or "minkowski".
#' @param hclust_method Cluster method to be used in `stats::hclust`. This
#'                      should be one of "ward.D", "ward.D2", "single", 
#'                      "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), 
#'                      "median" (= WPGMC) or "centroid" (= UPGMC).
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
Cluster <- function(
    Data,
    min_cluster_size,
    row_ids       = NULL,
    dist_method   = "euclidean",
    hclust_method = "ward.D2"
    ) {
  
  number_of_Samples  <- nrow(Data)
  number_of_Features <- ncol(Data) - 1
  
  # Warn if min_cluster_size is less than the number of features
  if (min_cluster_size < number_of_Features) {
    warning(sprintf(
      paste(
        "min_cluster_size (%d) is less than the number of ",
        "features/pigments (%d). This may lead to poor clustering or errors."
      ),
      min_cluster_size, number_of_Features
    ))
  }
  
  # Warn if min_cluster_size exceeds half the total number of samples
  maxAllowed <- floor(number_of_Samples / 2)
  if (min_cluster_size > maxAllowed) {
    warning(sprintf(
      paste(
        "min_cluster_size (%d) exceeds half of total samples (%d).", 
        "Clustering may not be meaningful."
      ),
      min_cluster_size, maxAllowed
    ))
  }
  
  # ---- safety & naming ---- #
  if (!is.null(row_ids)) {
    stopifnot(length(row_ids) == nrow(Data))
    rownames(Data) <- row_ids
  } else if (is.null(rownames(Data))) {
    rownames(Data) <- paste0("row_", seq_len(nrow(Data)))
  }

  # ---- helpers ---- #
  standardise <- function(Data) {
    b         <- Data[, -ncol(Data), drop = FALSE]
    b[b == 0] <- 1e-06
    b         <- b / Data[, ncol(Data)]
    v         <- sapply(b, bestNormalize::boxcox)
    v         <- do.call("cbind", v[1,])
    return(v)
  }

  ndf <- standardise(Data)
  rownames(ndf) <- rownames(Data)

  # ---- clustering ---- #
  mscluster <- stats::dist(ndf, method = dist_method)
  mv.hclust <- stats::hclust(mscluster, method = hclust_method)

  # make sure dendrogram shows your row_ids
  mv.hclust$labels <- rownames(ndf)

  # dynamic tree cut
  dynamicCut <- dynamicTreeCut::cutreeDynamic(
    mv.hclust,
    cutHeight            = 70,
    minClusterSize       = min_cluster_size,
    method               = "hybrid",
    distM                = as.matrix(stats::dist(ndf, method = dist_method)),
    deepSplit            = 4,
    pamStage             = TRUE,
    pamRespectsDendro    = TRUE,
    useMedoids           = FALSE,
    maxDistToLabel       = NULL,
    maxPamDist           = 50,
    respectSmallClusters = TRUE
  )

  # ---- outputs with row names preserved ----
  clust_list     <- split(Data, dynamicCut)

  assign_tbl <- data.frame(
    ID          = rownames(Data),
    Clust       = dynamicCut,
    row.names   = rownames(Data),
    check.names = FALSE
    )

  return(
    list(
      cluster.list = clust_list,
      cluster.plot = mv.hclust,
      assignments  = assign_tbl
      )
    )
}