#' Project data from one object into the UMAP embedding a single cell object
#' 
#' @description Adapted from: ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis
#' Jeffrey M. Granja, M. Ryan Corces, Sarah E. Pierce, S. Tansu Bagdatli, Hani Choudhry, Howard Y. Chang, William J. Greenleaf
#' doi: https://doi.org/10.1101/2020.04.28.066498
#' 
#' @param projector cell_data_set object with a reduced dimension matrix (currently iLSI supported) as specified in
#' reduced_dim argument and a model used to create a low dimensional embedding
#' @param projectee a SummarizedExperiment type object (cell_data_set currently supported) to be projected using 
#' the models contained in the projector
#' @param make_pseudo_single_cells whether to make pseudo-single cells from the data in the projectee (set this to true for bulk data)
#' @param ncells_coembedding number of cells in the projector to use in in the co-embedding with simulated single cells; default is 5000, will automatically
#' default to the total number of cells in the projector if less than this value.
#' @param reduced_dim A string specifying the reducedDim (currently LSI supported).
#' @param embedding A string specifying embedding type (currently UMAP supported).
#' @param n An integer specifying the number of subsampled "pseudo single cells" per bulk sample.  Note this is only relevant if 
#' make_pseudo_single_cells is TRUE
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to FALSE for a cleaner output.
#' @param threads The number of threads used for parallel execution
#' @export

project_data <- function(
  projector = NULL,
  projectee = NULL,
  ncells_coembedding = 5000,
  scale=F,
  reduced_dim = "LSI",
  embedding = "UMAP",
  make_pseudo_single_cells = FALSE,
  features = c("annotation-based", "range-based"),
  n = 250,
  verbose = TRUE,
  threads = 6,
  seed=2020
  
){
  check_input(input = projector, name = "projector", valid = c("cell_data_set"))
  check_input(input = projectee, name = "projectee", valid = c("SummarizedExperiment"))
  check_input(input = ncells_coembedding , name = "ncells_coembedding ", valid = c("numeric"))
  check_input(input = reduced_dim, name = "reduced_dim", valid = c("character"))
  check_input(input = embedding, name = "embedding", valid = c("character", "null"))
  check_input(input = n, name = "n", valid = c("integer"))
  check_input(input = verbose, name = "verbose", valid = c("boolean"))
  check_input(input = threads, name = "threads", valid = c("integer"))
  
  ##################################################
  # Extract data from bulk
  ##################################################
  rD<-reducedDims(projector)[[reduced_dim]]
  LSI_num_dim<-projector@preprocess_aux$iLSI$num_dim
  match(names(projector@reduce_dim_aux), embedding)
  embedding_num_dim<-projector@reduce_dim_aux[[embedding]]$num_dim
  
  sc_embedding<-reducedDims(projector)[[embedding]]
  rownames(sc_embedding)<-colnames(projector)
  
  if(features[1] %in% "annotation-based"){
    query<-projector@preprocess_aux$iLSI$features
    shared_rd<-extract_data(query, projectee)
  }
  
  if(features[1] %in% "range-based"){
    query<-projector@preprocess_aux$iLSI$granges
    shared_rd<-extract_data(query, projectee)
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
  if(make_pseudo_single_cells){
    depthN <- round(sum(projector@preprocess_aux$iLSI$row_sums / nrow(rD)))
    nRep <- 5
    n2 <- ceiling(n / nRep)
    ratios <- c(2, 1.5, 1, 0.5, 0.25) #range of ratios of number of fragments
  
    if(verbose) message(paste0("Simulating ", (n * dim(projectee)[2]), " single cells"))
    projRD <- pbmcapply::pbmclapply(seq_len(ncol(shared_rd$mat)), function(x){
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
      projRD <- as.matrix(projectLSI(simMat, LSI = projector@preprocess_aux$iLSI, verbose = verbose))
      rownames(projRD) <- paste0(colnames(shared_rd$mat)[x], "#", seq_len(nrow(projRD)))
      projRD
    }, mc.cores =  threads) %>% Reduce("rbind", .)
  
    # Deal with NaN
    if(any(is.nan(projRD))){
      projRD[is.nan(projRD)]<-0
      warning("NaN calculated during single cell generation")
    }
    
    if(scale){
      projRD <- scale_dims(projRD)
    }
  }else{
    projRD <- as.matrix(projectLSI(shared_rd$mat, LSI = projector@preprocess_aux$iLSI, verbose = verbose))
  }
  
  ##################################################
  # Check LSI and Embedding SVD columns
  ##################################################  
  
  
  if(embedding_num_dim != ncol(projRD)){
    stop("Error incosistency found with matching LSI dimensions to those used in embedding")
  }
  
  ##################################################
  # Get Previous UMAP Model
  ##################################################
  umap_model <- load_umap_model(projector@reduce_dim_aux[[embedding]]$model_file, embedding_num_dim)
  
  ##################################################
  # subsample
  ##################################################
  idx <- sort(sample(seq_len(nrow(rD)), min(nrow(rD), ncells_coembedding)))
  rD_ss <- rD[idx,,drop=FALSE]
  
  ##################################################
  # Project UMAP
  ##################################################
  if(verbose & make_pseudo_single_cells) message(paste0("Projecting simulated doublets onto manifold saved in: ", projector@reduce_dim_aux[[embedding]]$model_file))
  if(verbose & !make_pseudo_single_cells) message(paste0("Projecting projectee cells onto manifold saved in: ", projector@reduce_dim_aux[[embedding]]$model_file))
  set.seed(seed)
  #threads2 <- max(floor(threads/2), 1)
  simUMAP <- uwot::umap_transform(
    X = rbind(rD_ss, projRD), 
    model = umap_model, 
    verbose = verbose, 
    n_threads = threads
  )
  rownames(simUMAP) <- c(rownames(rD_ss), rownames(projRD))
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
  
  simUMAP <- DataFrame(simUMAP[rownames(projRD),,drop=FALSE])
  simUMAP$Type <- Rle(stringr::str_split(rownames(simUMAP), pattern = "#", simplify = TRUE)[,1])
  
  out <- SimpleList(
    simulatedBulkUMAP = simUMAP,
    singleCellUMAP = dfUMAP,
    simulatedReducedDims = projRD
  )
  return(out)
}