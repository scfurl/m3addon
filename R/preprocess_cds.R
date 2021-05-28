#' Preprocess a cds to prepare for trajectory inference
#'
#' @description Most analyses (including trajectory inference, and clustering)
#' in Monocle3, require various normalization and preprocessing steps.
#' \code{preprocess_cds} executes and stores these preprocessing steps.
#'
#' Specifically, depending on the options selected, \code{preprocess_cds} first
#' normalizes the data by log and size factor to address depth differences, or
#' by size factor only. Next, \code{preprocess_cds} calculates a lower
#' dimensional space that will be used as the input for further dimensionality
#' reduction like tSNE and UMAP.
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param method a string specifying the initial dimension method to use,
#'   currently either PCA or LSI. For LSI (latent semantic indexing), it
#'   converts the (sparse) expression matrix into tf-idf matrix and then
#'   performs SVD to decompose the gene expression / cells into certain
#'   modules / topics. Default is "PCA".
#' @param num_dim the dimensionality of the reduced space.
#' @param norm_method Determines how to transform expression values prior to
#'   reducing dimensionality. Options are "log", "size_only", and "none".
#'   Default is "log". Users should only use "none" if they are confident that
#'   their data is already normalized.
#' @param use_genes NULL or a list of gene IDs. If a list of gene IDs, only
#'   this subset of genes is used for dimensionality reduction. Default is
#'   NULL.
#' @param pseudo_count NULL or the amount to increase expression values before
#'   normalization and dimensionality reduction. If NULL (default), a
#'   pseudo_count of 1 is added for log normalization and 0 is added for size
#'   factor only normalization.
#' @param scaling When this argument is set to TRUE (default), it will scale
#'   each gene before running trajectory reconstruction. Relevant for
#'   method = PCA only.
#' @param verbose Whether to emit verbose output during dimensionality
#'   reduction
#' @param ... additional arguments to pass to limma::lmFit if
#'   residual_model_formula is not NULL
#' @return an updated cell_data_set object

preprocess_cds <- function(cds, method = c('PCA', "LSI"),
                           num_dim=50,
                           norm_method = c("log", "size_only", "none"),
                           use_genes = NULL,
                           pseudo_count=NULL,
                           scaling = TRUE,
                           verbose=FALSE,
                           ...) {
  
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be one of 'PCA' or 'LSI'")
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(norm_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "norm_method must be one of 'log', 'size_only' or 'none'")
  assertthat::assert_that(assertthat::is.count(num_dim))
  if(!is.null(use_genes)) {
    assertthat::assert_that(is.character(use_genes))
    assertthat::assert_that(all(use_genes %in% row.names(rowData(cds))),
                            msg = paste("use_genes must be NULL, or all must",
                                        "be present in the row.names of rowData(cds)"))
  }
  assertthat::assert_that(!is.null(size_factors(cds)),
                          msg = paste("You must call estimate_size_factors before calling",
                                      "preprocess_cds."))
  assertthat::assert_that(sum(is.na(size_factors(cds))) == 0,
                          msg = paste("One or more cells has a size factor of",
                                      "NA."))
  
  method <- match.arg(method)
  norm_method <- match.arg(norm_method)
  
  #ensure results from RNG sensitive algorithms are the same on all calls
  set.seed(2016)
  FM <- monocle3:::normalize_expr_data(cds, norm_method, pseudo_count)
  
  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }
  
  if (!is.null(use_genes)) {
    FM <- FM[use_genes, ]
  }else{
    use_genes <- rownames(FM)
  }
  
  fm_rowsums = Matrix::rowSums(FM)
  FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]
  
  if(method == 'PCA') {
    if (verbose) message("Remove noise by PCA ...")
    
    irlba_res <- m3addon:::sparse_prcomp_irlba(Matrix::t(FM),
                                     n = min(num_dim,min(dim(FM)) - 1),
                                     center = scaling, scale. = scaling)
    preproc_res <- irlba_res$x
    row.names(preproc_res) <- colnames(cds)
    
    irlba_rotation <- irlba_res$rotation
    row.names(irlba_rotation) <- rownames(FM)
    cds@preprocess_aux[[method]]$features <- use_genes
    cds@preprocess_aux[[method]]$svd<-list()
    cds@preprocess_aux[[method]]$svd$v <- irlba_res$rotation
    cds@preprocess_aux[[method]]$svd$d <- irlba_res$d
    cds@preprocess_aux[[method]]$svd$u <- irlba_res$u
    cds@preprocess_aux[[method]]$row_sums <- Matrix::rowSums(FM)
    cds@preprocess_aux[[method]]$gene_loadings <- irlba_rotation %*% diag(irlba_res$d)
    cds@preprocess_aux[[method]]$prop_var_expl <- irlba_res$sdev^2 / sum(irlba_res$sdev^2)
    svdDiag <- matrix(0, nrow=num_dim, ncol=num_dim)
    diag(svdDiag) <- cds@preprocess_aux[[method]]$svd$d
    matSVD <- t(svdDiag %*% t(irlba_res$rotation))
    rownames(matSVD) <- rownames(FM)
    colnames(matSVD) <- seq_len(ncol(matSVD))
    cds@preprocess_aux[[method]]$matSVD<-matSVD
    
  } else if(method == "LSI") {
    
    preproc_res <- tfidf(FM)
    num_col <- ncol(preproc_res)
    irlba_res <- irlba::irlba(Matrix::t(preproc_res),
                              nv = min(num_dim,min(dim(FM)) - 1))
    
    preproc_res <- irlba_res$u %*% diag(irlba_res$d)
    row.names(preproc_res) <- colnames(cds)
    
    irlba_rotation = irlba_res$v
    row.names(irlba_rotation) = rownames(FM)
    cds@preprocess_aux[[method]]$gene_loadings = irlba_rotation %*% diag(irlba_res$d)
    cds@preprocess_aux[[method]]$features <- use_genes
  }
  
  row.names(preproc_res) <- colnames(cds)
  
  reducedDims(cds)[[method]] <- as.matrix(preproc_res)
  cds@preprocess_aux[[method]]$beta = NULL
  
  cds
}


sparse_prcomp_irlba<-function (x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE, 
          ...) 
{
  a <- names(as.list(match.call()))
  ans <- list(scale = scale.)
  if ("tol" %in% a) 
    warning("The `tol` truncation argument from `prcomp` is not supported by\n            `prcomp_irlba`. If specified, `tol` is passed to the `irlba`\n            function to control that algorithm's convergence tolerance. See\n            `?prcomp_irlba` for help.")
  orig_x <- x
  if (class(x) != "DelayedMatrix") 
    x = DelayedArray::DelayedArray(x)
  args <- list(A = orig_x, nv = n)
  if (is.logical(center)) {
    if (center) 
      args$center <- DelayedMatrixStats::colMeans2(x)
  }
  else args$center <- center
  if (is.logical(scale.)) {
    if (is.numeric(args$center)) {
      scale. <- sqrt(DelayedMatrixStats::colVars(x))
      if (ans$scale) 
        ans$totalvar <- ncol(x)
      else ans$totalvar <- sum(scale.^2)
    }
    else {
      if (ans$scale) {
        scale. <- sqrt(DelayedMatrixStats::colSums2(x^2)/(max(1, 
                                                              nrow(x) - 1L)))
        ans$totalvar <- sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.)^2)/(nrow(x) - 
                                                                                   1L)))
      }
      else {
        ans$totalvar <- sum(DelayedMatrixStats::colSums2(x^2)/(nrow(x) - 
                                                                 1L))
      }
    }
    if (ans$scale) 
      args$scale <- scale.
  }
  else {
    args$scale <- scale.
    ans$totalvar <- sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.)^2)/(nrow(x) - 
                                                                               1L)))
  }
  if (!missing(...)) 
    args <- c(args, list(...))
  s <- do.call(irlba::irlba, args = args)
  ans$sdev <- s$d/sqrt(max(1, nrow(x) - 1))
  ans$rotation <- s$v
  ans$u<- s$u
  ans$d<- s$d
  colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), 
                                  sep = "")
  ans$center <- args$center
  if (retx) {
    ans <- c(ans, list(x = sweep(s$u, 2, s$d, FUN = `*`)))
    colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)), 
                             sep = "")
  }
  class(ans) <- c("irlba_prcomp", "prcomp")
  ans
}