#' Project Bulk RNA-seq data into single cell subspace
#' 
#' @description This function will Project Bulk RNA-seq data into single cell subspace. Adapted from: ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis
#' Jeffrey M. Granja, M. Ryan Corces, Sarah E. Pierce, S. Tansu Bagdatli, Hani Choudhry, Howard Y. Chang, William J. Greenleaf
#' doi: https://doi.org/10.1101/2020.04.28.066498
#' 
#' @param cds cell_data_set object
#' @param se a bulk Summarized Experiment.
#' @param reduced_dim A string specifying the reducedDim (currently LSI supported).
#' @param embedding A string specifying embedding type (UMAP supported).
#' @param n An integer specifying the number of subsampled "pseudo single cells" per bulk sample.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to FALSE for a cleaner output.
#' @param threads The number of threads used for parallel execution
#' @export

project_bulk_data <- function(
  cds = NULL,
  se = NULL,
  embedding_ncells = 5000,
  scale=F,
  reduced_dim = "LSI",
  embedding = "UMAP",
  features = c("annotation-based", "range-based"),
  n = 250,
  verbose = TRUE,
  threads = 6,
  force=F,
  binarize=F,
  seed=2020
  
){
  checkInput(input = cds, name = "cds", valid = c("cell_data_set"))
  checkInput(input = se, name = "se", valid = c("SummarizedExperiment"))
  checkInput(input = reduced_dim, name = "reduced_dim", valid = c("character"))
  checkInput(input = embedding, name = "embedding", valid = c("character", "null"))
  checkInput(input = n, name = "n", valid = c("integer"))
  checkInput(input = verbose, name = "verbose", valid = c("boolean"))
  checkInput(input = threads, name = "threads", valid = c("integer"))
  checkInput(input = binarize, name = "binarize", valid = c("boolean"))
  checkInput(input = force, name = "force", valid = c("boolean"))
  
  
  
  ##################################################
  # Reduced Dimensions
  ##################################################
  rD <- reducedDims(cds)[[reduced_dim]]
  if(features[1] %in% "annotation-based"){
    sub_se <- se[rownames(se) %in% cds@preprocess_aux$iLSI$features,]
    overlap <- dim(sub_se)[1]
    rows_in_bulk = cds@preprocess_aux$iLSI$features[cds@preprocess_aux$iLSI$features %in% row.names(sub_se)]
    bulk_mat<-as.matrix(getAssay(sub_se[rows_in_bulk,]))
    rownames(bulk_mat)<-rows_in_bulk
    total_sc_features=length(cds@preprocess_aux$iLSI$features)
    subset_rows<-cds@preprocess_aux$iLSI$features
  }
  #debug(getAssay)
  if(features[1] %in% "range-based"){
    rDGR <- cds@rowRanges[rownames(cds) %in% cds@preprocess_aux$iLSI$features]
    scfeat<-cds@preprocess_aux$iLSI$features
    # if("end" %in% colnames(rDFeatures)){
    #   rDGR <- GRanges(seqnames=rDFeatures$seqnames,IRanges(start=rDFeatures$start, end=rDFeatures$end))
    # }else{
    #   rDGR <- GRanges(seqnames=rDFeatures$seqnames,IRanges(start=rDFeatures$start, width = (rDFeatures$start) / (rDFeatures$idx - 1)))
    # }
    sub_se <- subsetByOverlaps(se, rDGR, ignore.strand = TRUE)
    #sub_se <- sub_se[order(rowSums(as.matrix(getAssay(sub_se, "counts"))), decreasing = TRUE), ]
    o <- DataFrame(findOverlaps(sub_se, rDGR, ignore.strand = TRUE))
    overlap <- length(unique(o[,2]))
    o <- o[!duplicated(o$subjectHits),]
    sub_se<-sub_se[o$queryHits, ]
    bulkfeat<-rownames(sub_se)
    rownames(sub_se) <- paste0("f", o$subjectHits)
    bulk_mat <- as.matrix(getAssay(sub_se, "counts"))
    rownames(bulk_mat)<-paste0("f", o$subjectHits)
    subset_rows = paste0("f", seq_along(rDGR))
    total_sc_features=length(rDGR)
    bulk_mat<-round(bulk_mat)
  }
  
  if(overlap == 0){
    stop(paste0("No overlaps between bulk RNA data and reduce dimensions feature found.",
                "\nEither recreate counts matrix or most likely these data sets are incompatible!"))
  }
  if( (overlap / total_sc_features) < 0.25 ){
    if(force){
      warning("Less than 25% of the features are present in this bulk RNA data set! Continuing since force = TRUE!")
    }else{
      stop("Less than 25% of the features are present in this bulk RNA data set! Set force = TRUE to continue!")
    }
  }
  
  ##################################################
  # Create Bulk Matrix
  ##################################################
  #undebug(getAssay)
  #undebug(safeSubset)
  bulkMat <- safeSubset(
    mat = bulk_mat, 
    subsetRows = subset_rows)
  #bm_sf<<-bulkMat

  message(paste0("Overlap Ratio of Reduced Dims Features = ", round(overlap / total_sc_features, 3)))
  #dim(bulkMat)
  cds@preprocess_aux$iLSI$features[1]
  exprs(cds["chr1_2082673_2083173",])
  ##################################################
  # Simulate and Project
  ##################################################
  depthN <- round(sum(cds@preprocess_aux$iLSI$row_sums / cds@preprocess_aux$iLSI$nCol))
  nRep <- 5
  n2 <- ceiling(n / nRep)
  ratios <- c(2, 1.5, 1, 0.5, 0.25) #range of ratios of number of fragments
  
  if(verbose) message(paste0("Simulating ", (n * dim(sub_se)[2]), " single cells"))
  simRD <- pbmcapply::pbmclapply(seq_len(ncol(bulkMat)), function(x){
    counts <- bulkMat[, x]
    counts <- rep(seq_along(counts), counts)
    simMat <- lapply(seq_len(nRep), function(y){
      ratio <- ratios[y]
      simMat <- matrix(sample(x = counts, size = ceiling(ratio * depthN) * n2, replace = TRUE), ncol = n2)
      simMat <- Matrix::summary(as(simMat, "dgCMatrix"))[,-1,drop=FALSE]
      simMat[,1] <- simMat[,1] + (y - 1) * n2
      simMat
    }) %>%  Reduce("rbind", .)
    simMat <- Matrix::sparseMatrix(i = simMat[,2], j = simMat[,1], x = rep(1, nrow(simMat)), dims = c(nrow(bulkMat), n2 * nRep))
    simRD <- as.matrix(projectLSI(simMat, LSI = cds@preprocess_aux$iLSI, verbose = verbose))
    rownames(simRD) <- paste0(colnames(bulkMat)[x], "#", seq_len(nrow(simRD)))
    simRD
  }, mc.cores =  threads) %>% Reduce("rbind", .)
  
  if(is.null(embedding)){
    if(rD$scaleDims){
      simRD <- scale_dims(simRD)
    }
    out <- SimpleList(
      simulatedReducedDims = simRD
    )
    return(out)
  }
  
  ##################################################
  # Prep Reduced Dims
  ##################################################
  sc_embedding<-list()
  sc_embedding$df <- cds@reduce_dim_aux[[embedding]]$embedding
  rownames(sc_embedding$df)<-colnames(cds)
  corCutOff <- 0.75
  #scaleDims <- cds@reduce_dim_aux[[embedding]]$scale_to
  dimsToUse <- cds@preprocess_aux$iLSI$num_dim #may need to fix this later...
  
  # if(is.null(scaleDims)){
  #   simRD <- scale_dims(simRD)
  # }
  
  if(scale){
    simRD <- scale_dims(simRD)
  }
  
  
  if(dimsToUse != ncol(simRD)){
    
    if(is.null(dimsToUse)){
      dimsToUse <- seq_len(ncol(rD))
    }
    
    # if(!is.null(corCutOff)){
    #   if(scaleDims){
    #     corToDepth <- rD$corToDepth$scaled
    #     dimsToUse <- dimsToUse[corToDepth < corCutOff]
    #   }else{
    #     corToDepth <- rD$corToDepth$none
    #     dimsToUse <- dimsToUse[corToDepth < corCutOff]
    #   }
    # }
    
    # if(embedding$params$nc != ncol(simRD)){
    #   if(verbose) message("Error incosistency found with matching LSI dimensions to those used in addEmbedding",
    #                       "\nReturning with simulated reduced dimension coordinates...")
    #   out <- SimpleList(
    #     simulatedReducedDims = simRD
    #   )
    #   return(out)
    # }
    
    simRD <- simRD[, dimsToUse, drop = FALSE]
    
  }
  #ArchR:::addUMAP
  #ArchR:::.saveUWOT
  #ArchR:::addUMAP
  ##################################################
  # Get Previous UMAP Model
  ##################################################
  umap_model <- load_umap_model(cds@reduce_dim_aux[[embedding]]$model_file, dimsToUse)
  
  idx <- sort(sample(seq_len(nrow(rD)), min(nrow(rD), embedding_ncells))) #Try to use 5000 or total cells to check validity
  rD_ss <- rD[idx,,drop=FALSE]
  
  ##################################################
  # Project UMAP
  ##################################################
  set.seed(1)
  #threads2 <- max(floor(threads/2), 1)
  simUMAP <- uwot::umap_transform(
    X = rbind(rD_ss, simRD), 
    model = umap_model, 
    verbose = verbose, 
    n_threads = threads
  )
  rownames(simUMAP) <- c(rownames(rD_ss), rownames(simRD))
  #Check if the projection matches using previous data
  c1 <- cor(simUMAP[rownames(rD_ss), 1], sc_embedding[[1]][rownames(rD_ss),1])
  c2 <- cor(simUMAP[rownames(rD_ss), 2], sc_embedding[[1]][rownames(rD_ss),2])
  if(min(c1, c2) < 0.8){
    message(paste0("Warning projection correlation is less than 0.8 (R = ", round(min(c1,c2), 4),").\nThese results may not be accurate because of the lack of heterogeneity in the single cell data."))
  }
  
  dfUMAP <- sc_embedding$df
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
checkInput<-function (input = NULL, name = NULL, valid = NULL) 
{
  valid <- unique(valid)
  if (is.character(valid)) {
    valid <- tolower(valid)
  }
  else {
    stop("Validator must be a character!")
  }
  if (!is.character(name)) {
    stop("name must be a character!")
  }
  if ("null" %in% tolower(valid)) {
    valid <- c("null", valid[which(tolower(valid) != "null")])
  }
  av <- FALSE
  for (i in seq_along(valid)) {
    vi <- valid[i]
    if (vi == "integer" | vi == "wholenumber") {
      if (all(is.numeric(input))) {
        cv <- min(abs(c(input%%1, input%%1 - 1))) < .Machine$double.eps^0.5
      }
      else {
        cv <- FALSE
      }
    }
    else if (vi == "null") {
      cv <- is.null(input)
    }
    else if (vi == "bool" | vi == "boolean" | vi == "logical") {
      cv <- is.logical(input)
    }
    else if (vi == "numeric") {
      cv <- is.numeric(input)
    }
    else if (vi == "vector") {
      cv <- is.vector(input)
    }
    else if (vi == "matrix") {
      cv <- is.matrix(input)
    }
    else if (vi == "sparsematrix") {
      cv <- is(input, "dgCMatrix")
    }
    else if (vi == "character") {
      cv <- is.character(input)
    }
    else if (vi == "factor") {
      cv <- is.factor(input)
    }
    else if (vi == "cell_data_set") {
      cv <- is(input, "cell_data_set")
    }
    else if (vi == "rlecharacter") {
      cv1 <- is(input, "Rle")
      if (cv1) {
        cv <- is(input@values, "factor") || is(input@values, 
                                               "character")
      }
      else {
        cv <- FALSE
      }
    }
    else if (vi == "palette") {
      cv <- all(.isColor(input))
    }
    else if (vi == "timestamp") {
      cv <- is(input, "POSIXct")
    }
    else if (vi == "dataframe" | vi == "data.frame" | vi == 
             "df") {
      cv1 <- is.data.frame(input)
      cv2 <- is(input, "DataFrame")
      cv <- any(cv1, cv2)
    }
    else if (vi == "fileexists") {
      cv <- all(file.exists(input))
    }
    else if (vi == "direxists") {
      cv <- all(dir.exists(input))
    }
    else if (vi == "granges" | vi == "gr") {
      cv <- is(input, "GRanges")
    }
    else if (vi == "grangeslist" | vi == "grlist") {
      cv <- .isGRList(input)
    }
    else if (vi == "list" | vi == "simplelist") {
      cv1 <- is.list(input)
      cv2 <- is(input, "SimpleList")
      cv <- any(cv1, cv2)
    }
    else if (vi == "bsgenome") {
      cv1 <- is(input, "BSgenome")
      cv2 <- tryCatch({
        library(input)
        eval(parse(text = input))
      }, error = function(e) {
        FALSE
      })
      cv <- any(cv1, cv2)
    }
    else if (vi == "se" | vi == "summarizedexperiment") {
      cv <- is(input, "SummarizedExperiment")
    }
    else if (vi == "seurat" | vi == "seuratobject") {
      cv <- is(input, "Seurat")
    }
    else if (vi == "txdb") {
      cv <- is(input, "TxDb")
    }
    else if (vi == "orgdb") {
      cv <- is(input, "OrgDb")
    }
    else if (vi == "bsgenome") {
      cv <- is(input, "BSgenome")
    }
    else if (vi == "parallelparam") {
      cv <- is(input, "BatchtoolsParam")
    }
    else {
      stop("Validator is not currently supported")
    }
    if (cv) {
      av <- TRUE
      break
    }
  }
  if (av) {
    return(invisible(TRUE))
  }
  else {
    stop("Input value for '", name, "' is not a ", paste(valid, 
                                                         collapse = ","), ", (", name, " = ", class(input), 
         ") please supply valid input!")
  }
}

#' @export
safeSubset<-function (mat = NULL, subsetRows = NULL, subsetCols = NULL) 
{
  if (!is.null(subsetRows)) {
    idxNotIn <- which(!subsetRows %in% rownames(mat))
    if (length(idxNotIn) > 0) {
      subsetNamesNotIn <- subsetRows[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, 
                                       dims = c(length(idxNotIn), ncol = ncol(mat)))
      dim(matNotIn)
      dim(mat)
      rownames(matNotIn) <- subsetNamesNotIn
      mat <- rbind(mat, matNotIn)
    }
    mat <- mat[subsetRows, ]
  }
  if (!is.null(subsetCols)) {
    idxNotIn <- which(subsetCols %ni% colnames(mat))
    if (length(idxNotIn) > 0) {
      subsetNamesNotIn <- subsetCols[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, 
                                       dims = c(nrow(mat), ncol = length(idxNotIn)))
      colnames(matNotIn) <- subsetNamesNotIn
      mat <- cbind(mat, matNotIn)
    }
    mat <- mat[, subsetCols]
  }
  mat
}

#' @export
getAssay<-function (se = NULL, assayName = NULL) 
{
  .assayNames <- function(se) {
    names(SummarizedExperiment::assays(se))
  }
  if (is.null(assayName)) {
    o <- SummarizedExperiment::assay(se)
  }
  else if (assayName %in% .assayNames(se)) {
    o <- SummarizedExperiment::assays(se)[[assayName]]
  }
  else {
    stop(sprintf("assayName '%s' is not in assayNames of se : %s", 
                 assayName, paste(.assayNames(se), collapse = ", ")))
  }
  return(o)
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
      idf   <- as(log(1 + LSI$nCol / LSI$row_sums), "sparseVector")
      if(verbose) message("Computing TF-IDF Matrix")
      mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*% 
        mat
    }
    else if (LSI$LSI_method == 2) {
      if(verbose) message("Computing Inverse Document Frequency")
      idf   <- as(LSI$nCol / LSI$row_sums, "sparseVector")
      if(verbose) message("Computing TF-IDF Matrix")
      mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*% 
        mat
      mat@x <- log(mat@x * LSI$scaleTo + 1)
    }else if (LSI$LSI_method == 3) {
      mat@x <- log(mat@x + 1)
      if(verbose) message("Computing Inverse Document Frequency")
      idf <- as(log(1 + LSI$nCol/LSI$row_sums), "sparseVector")
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
