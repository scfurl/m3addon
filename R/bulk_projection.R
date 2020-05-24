#' Project Bulk RNA-seq data into single cell subspace
#' 
#' @description This function will Project bulk sequencing data into single cell subspace. Adapted from: ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis
#' Jeffrey M. Granja, M. Ryan Corces, Sarah E. Pierce, S. Tansu Bagdatli, Hani Choudhry, Howard Y. Chang, William J. Greenleaf
#' doi: https://doi.org/10.1101/2020.04.28.066498
#' 
#' @param cds cell_data_set object
#' @param se a bulk Summarized Experiment.
#' @param ncells_coembedding number of true single cells to use in in the co-embedding with simulated single cells 
#' @param reduced_dim A string specifying the reducedDim (currently LSI supported).
#' @param embedding A string specifying embedding type (currently UMAP supported).
#' @param n An integer specifying the number of subsampled "pseudo single cells" per bulk sample.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to FALSE for a cleaner output.
#' @param threads The number of threads used for parallel execution
#' @export

project_bulk_data <- function(
  cds = NULL,
  se = NULL,
  ncells_coembedding = 5000,
  scale=F,
  reduced_dim = "LSI",
  embedding = "UMAP",
  features = c("annotation-based", "range-based"),
  n = 250,
  verbose = TRUE,
  threads = 6,
  seed=2020
  
){
  check_input(input = cds, name = "cds", valid = c("cell_data_set"))
  check_input(input = se, name = "se", valid = c("SummarizedExperiment"))
  check_input(input = ncells_coembedding , name = "ncells_coembedding ", valid = c("numeric"))
  check_input(input = reduced_dim, name = "reduced_dim", valid = c("character"))
  check_input(input = embedding, name = "embedding", valid = c("character", "null"))
  check_input(input = n, name = "n", valid = c("integer"))
  check_input(input = verbose, name = "verbose", valid = c("boolean"))
  check_input(input = threads, name = "threads", valid = c("integer"))
  
  ##################################################
  # Extract data from bulk
  ##################################################
  rD<-reducedDims(cds)[[reduced_dim]]
  LSI_num_dim<-cds@preprocess_aux$iLSI$num_dim
  match(names(cds@reduce_dim_aux), embedding)
  embedding_num_dim<-cds@reduce_dim_aux[[embedding]]$num_dim
  
  sc_embedding<-reducedDims(cds)[[embedding]]
  rownames(sc_embedding)<-colnames(cds)

  if(features[1] %in% "annotation-based"){
    query<-cds@preprocess_aux$iLSI$features
    shared_rd<-extract_data(query, se)
  }
  
  if(features[1] %in% "range-based"){
    query<-cds@preprocess_aux$iLSI$Granges
    shared_rd<-extract_data(query, se)
  }
  
  message(paste0("Overlap Ratio of Reduced Dims Features = ", round(shared_rd$overlap, 3)))
  
  if( (shared_rd$overlap) < 0.25 ){
    if(force){
      warning("Less than 25% of the features are present in this bulk RNA data set! Continuing since force = TRUE!")
    }else{
      stop("Less than 25% of the features are present in this bulk RNA data set! Set force = TRUE to continue!")
    }
  }
  
  ##################################################
  # Simulate single cells and project using original LSI/SVD model
  ##################################################
  depthN <- round(sum(cds@preprocess_aux$iLSI$row_sums / nrow(rD)))
  nRep <- 5
  n2 <- ceiling(n / nRep)
  ratios <- c(2, 1.5, 1, 0.5, 0.25) #range of ratios of number of fragments
  
  if(verbose) message(paste0("Simulating ", (n * dim(se)[2]), " single cells"))
  simRD <- pbmcapply::pbmclapply(seq_len(ncol(shared_rd$mat)), function(x){
    counts <- shared_rd$mat[, x]
    counts <- rep(seq_along(counts), counts)
    simMat <- lapply(seq_len(nRep), function(y){
      ratio <- ratios[y]
      simMat <- matrix(sample(x = counts, size = ceiling(ratio * depthN) * n2, replace = TRUE), ncol = n2)
      simMat <- Matrix::summary(as(simMat, "dgCMatrix"))[,-1,drop=FALSE]
      simMat[,1] <- simMat[,1] + (y - 1) * n2
      simMat
    }) %>%  Reduce("rbind", .)
    simMat <- Matrix::sparseMatrix(i = simMat[,2], j = simMat[,1], x = rep(1, nrow(simMat)), dims = c(nrow(shared_rd$mat), n2 * nRep))
    simRD <- as.matrix(projectLSI(simMat, LSI = cds@preprocess_aux$iLSI, verbose = verbose))
    rownames(simRD) <- paste0(colnames(shared_rd$mat)[x], "#", seq_len(nrow(simRD)))
    simRD
  }, mc.cores =  threads) %>% Reduce("rbind", .)
  
  # Deal with NaN
  if(any(is.nan(simRD))){
    simRD[is.nan(simRD)]<-0
    warning("NaN calculated during single cell generation")
  }
  
  if(scale){
    simRD <- scale_dims(simRD)
  }
  
  
  ##################################################
  # Check LSI and Embedding SVD columns
  ##################################################


  if(embedding_num_dim != ncol(simRD)){
    stop("Error incosistency found with matching LSI dimensions to those used in embedding")
  }
  
  ##################################################
  # Get Previous UMAP Model
  ##################################################
  umap_model <- load_umap_model(cds@reduce_dim_aux[[embedding]]$model_file, embedding_num_dim)
  
  ##################################################
  # subsample
  ##################################################
  idx <- sort(sample(seq_len(nrow(rD)), min(nrow(rD), ncells_coembedding)))
  rD_ss <- rD[idx,,drop=FALSE]
  
  ##################################################
  # Project UMAP
  ##################################################
  if(verbose) message(paste0("Projecting simulated doublets onto manifold saved in: ", cds@reduce_dim_aux[[embedding]]$model_file))
  set.seed(seed)
  #threads2 <- max(floor(threads/2), 1)
  simUMAP <- uwot::umap_transform(
    X = rbind(rD_ss, simRD), 
    model = umap_model, 
    verbose = verbose, 
    n_threads = threads
  )
  rownames(simUMAP) <- c(rownames(rD_ss), rownames(simRD))
  ##################################################
  # Check correlation of subsampled cells
  ##################################################
  c1 <- cor(simUMAP[rownames(rD_ss), 1], sc_embedding[rownames(rD_ss),1])
  c2 <- cor(simUMAP[rownames(rD_ss), 2], sc_embedding[rownames(rD_ss),2])
  if(min(c1, c2) < 0.8){
    message(paste0("Warning projection correlation is less than 0.8 (R = ", round(min(c1,c2), 4),").\nThese results may not be accurate because of the lack of heterogeneity in the single cell data."))
  }
  
  dfUMAP <- sc_embedding
  colnames(dfUMAP) <- c("UMAP1", "UMAP2")
  colnames(simUMAP) <- c("UMAP1", "UMAP2")
  dfUMAP <- DataFrame(dfUMAP)
  dfUMAP$Type <- Rle("single_cell", lengths = nrow(dfUMAP))
  
  simUMAP <- DataFrame(simUMAP[rownames(simRD),,drop=FALSE])
  simUMAP$Type <- Rle(stringr::str_split(rownames(simUMAP), pattern = "#", simplify = TRUE)[,1])
  
  out <- SimpleList(
    simulatedBulkUMAP = simUMAP,
    singleCellUMAP = dfUMAP,
    simulatedReducedDims = simRD
  )
  return(out)
}


#' @export
projectLSI<-function (mat = NULL, LSI = NULL, returnModel = FALSE, verbose = FALSE, seed=2020) 
{
  out2 <- tryCatch({
    require(Matrix)
    set.seed(seed)
    if(verbose) message(sprintf("Projecting LSI, Input Matrix = %s GB", 
                         round(object.size(mat)/10^9, 3)))
    if(verbose) message("Subsetting by Non-Zero features in inital Matrix")
    #mat <- mat[LSI$idx, ] I don't see what this does after safeSubset
    if (LSI$binarize) {
      if(verbose) message("Binarizing Matrix")
      mat@x[mat@x > 0] <- 1
    }
    if(verbose) message("Computing Term Frequency")
    colSm <- Matrix::colSums(mat)
    if (any(colSm == 0)) {
      exclude <- which(colSm == 0)
      mat <- mat[, -exclude]
      colSm <- colSm[-exclude]
    }
    mat@x <- mat@x/rep.int(colSm, Matrix::diff(mat@p))
    if (LSI$LSI_method == 1) {
      if(verbose) message("Computing Inverse Document Frequency")
      idf   <- as(log(1 + nrow(LSI$svd$v) / LSI$row_sums), "sparseVector")
      if(verbose) message("Computing TF-IDF Matrix")
      mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*% 
        mat
    }
    else if (LSI$LSI_method == 2) {
      if(verbose) message("Computing Inverse Document Frequency")
      idf   <- as( nrow(LSI$svd$v) / LSI$row_sums, "sparseVector")
      if(verbose) message("Computing TF-IDF Matrix")
      mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*% 
        mat
      mat@x <- log(mat@x * LSI$scaleTo + 1)
    }else if (LSI$LSI_method == 3) {
      mat@x <- log(mat@x + 1)
      if(verbose) message("Computing Inverse Document Frequency")
      idf <- as(log(1 + nrow(LSI$svd$v) /LSI$row_sums), "sparseVector")
      if(verbose) message("Computing TF-IDF Matrix")
      mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*% 
        mat
      }else {
      stop("LSIMethod unrecognized please select valid method!")
    }
    gc()
    idxNA <- Matrix::which(is.na(mat), arr.ind = TRUE)
    if (length(idxNA) > 0) {
      if(verbose) message((sprintf("Zeroing %s NA elements", length(idxNA))))
      mat[idxNA] <- 0
    }
    if(verbose) message("Calculating V Matrix")
    V <- Matrix::t(mat) %*% LSI$svd$u %*% Matrix::diag(1/LSI$svd$d)
    if(verbose) message("Computing Projected Coordinates")
    svdDiag <- matrix(0, nrow = LSI$num_dim, ncol = LSI$num_dim)
    diag(svdDiag) <- LSI$svd$d
    matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
    matSVD <- as.matrix(matSVD)
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI", seq_len(ncol(matSVD)))
    if (returnModel) {
      if(verbose) message("Calculating Re-Projected Matrix")
      X <- LSI$svd$u %*% diag(LSI$svd$d) %*% t(V)
      out <- list(matSVD = matSVD, V = V, X = X)
    }
    else {
      out <- matSVD
    }
    out
  }, error = function(e) {
    errorList <- list(mat = mat, colSm = if (exists("colSm", 
                                                    inherits = FALSE)) colSm else "Error with colSm!", 
                      idf = if (exists("idf", inherits = FALSE)) idf else "Error with idf!", 
                      V = if (exists("V", inherits = FALSE)) V else "Error with V!", 
                      matSVD = if (exists("matSVD", inherits = FALSE)) matSVD else "Error with matSVD!")
    errorLog(e, fn = "projectLSI", info = "", errorList = errorList)
  })
  out2
}

#' @export
scale_dims<-function(x, scale_max = NULL){
  if(!is.null(scale_max)){
    row_z_scores(m=x, min=-scale_max, max = scale_max, limit = TRUE)
  }else{
    row_z_scores(m=x)
  }
}

#' @export
row_z_scores<-function(m = NULL, min = -2, max = 2, limit = FALSE){
  z <- sweep(m - Matrix::rowMeans(m), 1, matrixStats::rowSds(m),`/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}

#' @export
errorLog<-function(
  e = NULL,
  fn = NULL,
  info = NULL, 
  errorList = NULL,
  throwError = TRUE
){
  
  header <- "************************************************************"
  if(!is.null(errorList)){
    tryCatch({
      saveRDS(errorList, "Save-Error.rds")
      message("Saving a list of errors to Save-Error.rds")
    }, error = function(e){
      message("Error recording errorList")
    })
  }
  print(e)
  cat(sprintf("\n%s\n\n", header))
  if(throwError) stop("Exiting See Error Above")
}
