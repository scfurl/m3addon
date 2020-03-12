#' Compute a projection of a cell_data_set object into a lower dimensional
#' space with non-linear dimension reduction methods
#'
#' @description Monocle3 aims to learn how cells transition through a
#' biological program of gene expression changes in an experiment. Each cell
#' can be viewed as a point in a high-dimensional space, where each dimension
#' describes the expression of a different gene. Identifying the program of
#' gene expression changes is equivalent to learning a \emph{trajectory} that
#' the cells follow through this space. However, the more dimensions there are
#' in the analysis, the harder the trajectory is to learn. Fortunately, many
#' genes typically co-vary with one another, and so the dimensionality of the
#' data can be reduced with a wide variety of different algorithms. Monocle3
#' provides two different algorithms for dimensionality reduction via
#' \code{reduce_dimensions} (UMAP and tSNE). The function
#' \code{reduce_dimensions} is the second step in the trajectory building
#' process after \code{preprocess_cds}.
#'
#' UMAP is implemented from the package uwot.
#'
#' @param cds the cell_data_set upon which to perform this operation.
#' @param max_components the dimensionality of the reduced space. Default is 2.
#' @param reduction_method A character string specifying the algorithm to use
#'   for dimensionality reduction. Currently "UMAP", "tSNE", and "PCA" are
#'   supported.
#'@param num_dim Numeric indicating the number of prinicipal components to be 
#'    in downstream ordering.  Default value is NULL which will result in use 
#'    of all PCs
#' @param preprocess_method A string indicating the preprocessing method used
#'   on the data. Options are "PCA" and "LSI". Default is "LSI".
#' @param umap.metric A string indicating the distance metric to be used when
#'   calculating UMAP. Default is "cosine". See uwot package's
#'   \code{\link[umap]{umap}} for details.
#' @param umap.min_dist Numeric indicating the minimum distance to be passed to
#'   UMAP function. Default is 0.1.See uwot package's \code{\link[umap]{umap}}
#'   for details.
#' @param umap.n_neighbors Integer indicating the number of neighbors to use
#'   during kNN graph construction. Default is 15L. See uwot package's
#'   \code{\link[umap]{umap}} for details.
#' @param umap.fast_sgd Logical indicating whether to use fast SGD. Default is
#'   TRUE. See uwot package's \code{\link[umap]{umap}} for details.
#' @param umap.nn_method String indicating the nearest neighbor method to be
#'   used by UMAP. Default is "annoy". See uwot package's
#'   \code{\link[umap]{umap}} for details.
#' @param cores Number of compute cores to use.
#' @param verbose Logical, whether to emit verbose output.
#' @param ... additional arguments to pass to the dimensionality reduction
#'   function.
#' @return an updated cell_data_set object
#' @references UMAP: McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation
#'   and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
#' @references tSNE: Laurens van der Maaten and Geoffrey Hinton. Visualizing
#'   data using t-SNE. J. Mach. Learn. Res., 9(Nov):2579– 2605, 2008.
#' @references monocle3 - this function only differs from that found in monocle3 in \
#' that it allows for selection of num_dim
#' @export
reduce_dimension <- function(cds,
                             max_components=2,
                             reduction_method=c("UMAP", 'tSNE', 'PCA'),
                             preprocess_method=c("PCA", "LSI"),
                             umap.metric = "cosine",
                             umap.min_dist = 0.1,
                             umap.n_neighbors = 15L,
                             umap.fast_sgd = FALSE,
                             umap.nn_method = "annoy",
                             cores=1,
                             verbose=FALSE,
                             num_dim=NULL,
                             ...){
  extra_arguments <- list(...)
  
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'UMAP', 'PCA' or 'tSNE'")
  
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(preprocess_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "preprocess_method must be one of 'PCA' or 'LSI'")
  
  reduction_method <- match.arg(reduction_method)
  preprocess_method <- match.arg(preprocess_method)
  
  assertthat::assert_that(assertthat::is.count(max_components))
  
  assertthat::assert_that(!is.null(reducedDims(cds)[[preprocess_method]]),
                          msg = paste("Data has not been preprocessed with",
                                      "chosen method:", preprocess_method,
                                      "Please run preprocess_cds with",
                                      "method =", preprocess_method,
                                      "before running reduce_dimension."))
  if(reduction_method == "PCA") {
    assertthat::assert_that(preprocess_method == "PCA",
                            msg = paste("preprocess_method must be 'PCA' when",
                                        "reduction_method = 'PCA'"))
    assertthat::assert_that(!is.null(reducedDims(cds)[["PCA"]]),
                            msg = paste("When reduction_method = 'PCA', the",
                                        "cds must have been preprocessed for",
                                        "PCA. Please run preprocess_cds with",
                                        "method = 'PCA' before running",
                                        "reduce_dimension with",
                                        "reduction_method = 'PCA'."))
  }
  
  #ensure results from RNG sensitive algorithms are the same on all calls
  set.seed(2016)
  
  if (reduction_method=="UMAP" && (umap.fast_sgd == TRUE || cores > 1)){
    message(paste("Note: reduce_dimension will produce slightly different",
                  "output each time you run it unless you set",
                  "'umap.fast_sgd = FALSE' and 'cores = 1'"))
  }
  
  preprocess_mat <- reducedDims(cds)[[preprocess_method]]
  if(!is.null(num_dim)){
    preprocess_mat<-preprocess_mat[,1:num_dim]
  }
  
  if(reduction_method == "PCA") {
    if (verbose) message("Returning preprocessed PCA matrix")
  } else if (reduction_method == "tSNE") {
    if (verbose) message("Reduce dimension by tSNE ...")
    
    tsne_res <- Rtsne::Rtsne(as.matrix(preprocess_mat), dims = max_components,
                             pca = F, check_duplicates=FALSE, ...)
    
    tsne_data <- tsne_res$Y[, 1:max_components]
    row.names(tsne_data) <- colnames(tsne_data)
    
    reducedDims(cds)$tSNE <- tsne_data
    
  } else if (reduction_method == c("UMAP")) {
    if (verbose)
      message("Running Uniform Manifold Approximation and Projection")
    
    umap_res = uwot::umap(as.matrix(preprocess_mat),
                          n_components = max_components,
                          metric = umap.metric,
                          min_dist = umap.min_dist,
                          n_neighbors = umap.n_neighbors,
                          fast_sgd = umap.fast_sgd,
                          n_threads=cores,
                          verbose=verbose,
                          nn_method = umap.nn_method,
                          ...)
    
    row.names(umap_res) <- colnames(cds)
    reducedDims(cds)$UMAP <- umap_res
  }
  
  ## Clear out any old graphs:
  cds@principal_graph_aux[[reduction_method]] <- NULL
  cds@principal_graph[[reduction_method]] <- NULL
  cds@clusters[[reduction_method]] <- NULL
  
  cds
}




#' Iterative LSI
#'
#' @description This function aims to both minimize batch effects and accentuate
#' cell type differences in a single cell experiment.  This function was implemented
#' using Monocle3 but takes inspiration from the Granja et. al. reference cited below which took inspiration from the fly ATAC paper. At it's heart
#' this function iterates through three main steps: 1) Using TFIDF transformation and SVD
#' to normalize data 2) Clustering this normalized data using leiden clustering in high dimensional
#' space and 3) identifying those features that are over-represented in the resulting clusters using
#' a simple counting method.  These three steps are repeated using features identified in step 3 to subset 
#' the normalization matrix in step 1 and repeating through the process.  TFIDF transformation is 
#' supplied in this package.  SVD is performed using the irilba package.  Leiden clustering is performed using
#' the monocle3 implementation and finally the counting per cluster is performed using the edgeR cpm function.  This
#' function takes as its input a cell_data_set and will iterate through n number of iterations.  The output of this function
#' is then appropriately input into dimensionality reduction methods such as UMAP or tSNE.  The number of iterations
#' is set by the number of resolution parameters specified.
#' 
#' @param cds the cell_data_set upon which to perform this operation.
#' @param num_dim Numeric indicating the number of prinicipal components to be 
#'    in downstream ordering.  Default value is NULL which will result in use 
#'    of all PCs
#' @param resolution vector of resolution values for leiden clustering
#' @param binarize boolean whether to binarize data prior to TFIDF transformation
#' @param nFeatures number of features to use for dimensionality reduction (default 3000).  To use different numbers
#' of features for different iterations, supply a vector that is the same length as the resolution vector.
#' @return an updated cell_data_set object with a reduced dimension LSI object and clusters object
#' @references Granja, J. M.et al. (2019). Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype 
#' acute leukemia. Nature Biotechnology, 37(12), 1458–1465.
#' @references UMAP: McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation
#'   and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
#' @references tSNE: Laurens van der Maaten and Geoffrey Hinton. Visualizing
#'   data using t-SNE. J. Mach. Learn. Res., 9(Nov):2579– 2605, 2008.
#' @references Cusanovich, D. A., Reddington, J. P., Garfield, D. A., Daza, R. M., Aghamirzaie, D., Marco-Ferreres, R., et al. (2018). The 
#'   cis-regulatory dynamics of embryonic development at single-cell resolution. Nature, 555(7697), 538–542.
#' @export
iterative_LSI <- function(cds,
                             num_dim=25,
                          resolution=c(1e-4, 3e-4, 5e-4),
                          nFeatures=c(3000,3000,3000), 
                          binarize=FALSE,
                          seed=2020, scaleTo = 10000, leiden_k=20, leiden_weight=FALSE, leiden_iter=1, verbose=F,
                          ...){
  extra_arguments <- list(...)
  if(length(resolution)<2){
    stop("This method is intended to iterate.  Adjust number of resolution elements and retry")
  }
  if(length(nFeatures)!=length(resolution)){
    message("Numbers of elements for resolution and nFeatures do not match.  Will use nFeatures[1]...")
    nFeatures<-rep(nFeatures, length(resolution))
  }
  mat<-assay(cds)
  set.seed(seed)
  if(binarize){
    message("Binarizing...")
    mat@x[mat@x > 0] <- 1
  }
  matNorm <- t(t(mat)/Matrix::colSums(mat)) * scaleTo
  matNorm@x <- log2(matNorm@x + 1)
  message("Performing LSI/SDF for iteration 1....")
  tf<-tf_idf_transform(mat[head(order(sparseRowVariances(matNorm),decreasing=TRUE), topN),])
  tf@x[is.na(tf@x)] <- 0
  matSVD<-svd_lsi(tf, num_dim)
  cluster_result <- monocle3:::leiden_clustering(data = matSVD, 
                                      pd = pData(cds), k = leiden_k, weight = leiden_weight, num_iter = leiden_iter, 
                                      resolution_parameter = resolution[1], random_seed = seed, 
                                      verbose = verbose, ...)
  clusterMat <- edgeR::cpm(groupSums(mat, factor(cluster_result$optim_res$membership), sparse = TRUE), log=TRUE, prior.count = 3)
  for(iterations in 2:length(resolution)){
    message("Performing LSI/SDF for iteration ", iterations, "....")
    tf<-tf_idf_transform(mat[head(order(rowVars(clusterMat), decreasing=TRUE), topN),])
    tf@x[is.na(tf@x)] <- 0
    if(iterations!=length(resolution)){
      matSVD<-svd_lsi(tf, num_dim, mat_only=TRUE)
      cluster_result <- monocle3:::leiden_clustering(data = matSVD, 
                                                     pd = pData(cds), k = leiden_k, weight = leiden_weight, num_iter = leiden_iter, 
                                                     resolution_parameter = resolution[iterations], random_seed = seed, 
                                                     verbose = verbose, ...)
      clusterMat <- edgeR::cpm(groupSums(mat, factor(cluster_result$optim_res$membership), sparse = TRUE), log=TRUE, prior.count = 3)
    }else{
      svd_list<-svd_lsi(tf, num_dim, mat_only=FALSE)
      reducedDims(cds)[["LSI"]]<-svd_list$matSVD
      irlba_rotation = svd_list$svd$v
      row.names(irlba_rotation) = colnames(cds)
      cds@preprocess_aux$gene_loadings = irlba_rotation
      cds@clusters[["LSI"]]<- monocle3:::leiden_clustering(data = svd_list$matSVD, 
                                                     pd = pData(cds), k = leiden_k, weight = leiden_weight, num_iter = leiden_iter, 
                                                     resolution_parameter = resolution[iterations], random_seed = seed, 
                                                     verbose = verbose, ...)
    }
  }
  cds
}
