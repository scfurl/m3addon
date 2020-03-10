#'Cluster cells using Density Peak algorithm
#'
#' @description
#' Unsupervised clustering of cells is a common step in many single-cell expression workflows. \
#' In an experiment containing a mixture of cell types, each cluster might correspond to a different \
#' cell type. This method takes a cell_data_set as input along with a requested number of clusters, \
#' clusters them using density peak clustering), and then returns the cell_data_set with the cluster assignments \
#' stored in the pData table as "Cluster". Use the plot_rho_delta to visualize the rho and delta parameters that will \
#' help determine the number of clusters to cluster.
#' @param cds the cell_data_set upon which to perform this operation
#' @param rho The threshold of local density (rho) used to select the density peaks
#' @param delta The threshold of local distance (delta) used to select the density peaks
#' @param gaussian A logic flag passed to densityClust function in desnityClust package to determine whether or not Gaussian kernel will be used for calculating the local density
#' @return an updated cell_data_set object
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. \
#' Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' @importFrom densityClust densityClust
#' @importFrom densityClust findClusters
#' @export
#' 
density_peak<-function(cds, rho=NULL, delta=NULL, reduction_method=c("UMAP", 'tSNE'), gaussian=T, pData_col="Cluster"){
  set.seed(2019)
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'UMAP' or 'tSNE'")
  cds@reduce_dim_aux$densityPeak<- densityClust(dist(reducedDims(cds)[[reduction_method]]), gaussian = gaussian)
  if (is.null(rho) | is.null(delta)) {
    message("Use 0.95 of the delta and 0.95 of the rho as the cutoff for assigning density peaks and clusters")
    rho <- quantile(cds@reduce_dim_aux$densityPeak$rho, probs = 0.95)
    delta <- quantile(cds@reduce_dim_aux$densityPeak$delta, probs = 0.95)
  }
  #plot_rho_delta(cds)
  dc<-findClusters(cds@reduce_dim_aux$densityPeak, rho=rho, delta=delta)
  pData(cds)[[pData_col]]<-factor(dc$clusters)
  cds
}
