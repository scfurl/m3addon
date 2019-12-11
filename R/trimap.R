#' Dimensionality Reduction Using Triplet Constraints
#' @description Find a low-dimensional repersentation of the data by satisfying the sampled triplet constraints from the high-dimensional features.
#' @param data input data samples (rows), features (columns)
#' @param n_dims Number of dimensions of the embedding (default = 2)
#' @param n_inliers: Number of inlier points for triplet constraints (default = 10)
#' @param n_outliers: Number of outlier points for triplet constraints (default = 5)
#' @param n_random: Number of random triplet constraints per point (default = 5)
#' @param distance: Distance measure ('euclidean' (default), 'manhattan', 'angular', 'hamming')
#' @param lr: Learning rate (default = 1000.0)
#' @param n_iters: Number of iterations (default = 400)
#' @param knn_tuple: Use the pre-computed nearest-neighbors information in form of a tuple (knn_nbrs, knn_distances) (default = None)
#' @param apply_pca: Apply PCA to reduce the dimensions to 100 if necessary before the nearest-neighbor calculation (default = True)
#' @param opt_method: Optimization method ('sd': steepest descent,  'momentum': GD with momentum, 
#' 'dbd': GD with momentum delta-bar-delta (default))
#' @param verbose: Print the progress report (default = True)
#' @param weight_adj: Adjusting the weights using a non-linear transformation (default = 500.0)
#' @param return_seq: Return the sequence of maps recorded every 10 iterations (default = False)
#' @export
#' @importFrom reticulate source_python
#' @references TriMap: Large-scale Dimensionality Reduction Using Triplets. E Amid, MK Warmuth - arXiv preprint arXiv:1910.00204, 2019 - arxiv.org
trimap<-function(cds, python_home = system("which python", intern = TRUE), num_dims = NULL, module_file=paste(system.file(package = "m3addon"),"trimap.py", sep = "/"),
                  n_dims = 2, n_inliers= 10, n_outliers = 5, n_random = 5, distance = c('euclidean', 'manhattan', 'angular', 'hamming'), 
                  lr = 1000.0, n_iters = 400, knn_tuple = NULL, apply_pca_trimap = FALSE, opt_method = c('dbd', 'sd',  'momentum'), 
                  verbose = TRUE, weight_adj = 500.0,  return_seq = FALSE)
{
  #reticulate::use_python(python_home)
  if(!py_available("trimap")) stop("python module trimap does not seem to be installed; - try running 'py_config()'")
  source_python(module_file)
  if(is.null(reducedDims(cds)["PCA"]))stop("PCA not found in the reducedDims slot")
  if(is.null(num_dims)){
    num_dims<-dim(reducedDims(cds)["PCA"][[1]])[2]
  }
  data <- reducedDims(cds)["PCA"][[1]][,1:num_dims]
  #data<-t(as.matrix(exprs(cds)))
  distance <- match.arg(distance)
  opt_method <- match.arg(opt_method)
  trimap_args<-c(list(data, n_dims, n_inliers, n_outliers, n_random, distance, lr, n_iters, knn_tuple, apply_pca_trimap, opt_method,
                                      verbose, weight_adj, return_seq))
  reducedDims(cds)[["trimap"]] <- do.call(trimap_fromR, trimap_args)
  return(cds)
}

