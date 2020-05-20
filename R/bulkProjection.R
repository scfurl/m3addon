#' Project Bulk RNA-seq data into single cell subspace
#' 
#' This function will Project Bulk RNA-seq data into single cell subspace. Adapted from: ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis
#' Jeffrey M. Granja, M. Ryan Corces, Sarah E. Pierce, S. Tansu Bagdatli, Hani Choudhry, Howard Y. Chang, William J. Greenleaf
#' doi: https://doi.org/10.1101/2020.04.28.066498
#' 
#' @param cds cell_data_set object
#' @param seRNA Bulk RNA Summarized Experiment.
#' @param reducedDims A string specifying the reducedDims.
#' @param embedding A string specifying embedding.
#' @param n An integer specifying the number of subsampled "pseudo single cells" per bulk sample.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to FALSE for a cleaner output.
#' @param threads The number of threads used for parallel execution
#' @export
#'
projectBulk <- function(
  cds = NULL,
  seRNA = NULL,
  reducedDims = "LSI",
  embedding = "UMAP",
  n = 250,
  verbose = TRUE,
  threads = 6,
  force=F,
  binarize=F,
  LSIMethod=1
){
  
  checkInput(input = cds, name = "cds", valid = c("cell_data_set"))
  checkInput(input = seRNA, name = "seRNA", valid = c("SummarizedExperiment"))
  checkInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  checkInput(input = embedding, name = "embedding", valid = c("character", "null"))
  checkInput(input = n, name = "n", valid = c("integer"))
  checkInput(input = verbose, name = "verbose", valid = c("boolean"))
  checkInput(input = threads, name = "threads", valid = c("integer"))
  checkInput(input = binarize, name = "binarize", valid = c("boolean"))
  checkInput(input = force, name = "force", valid = c("boolean"))
  checkInput(input = LSIMethod, name = "LSIMethod", valid = c("integer"))
  
  ##################################################
  # Reduced Dimensions
  ##################################################
  rD <- reducedDims(cds)[[reducedDims]]
  subRNA <- seRNA[rownames(seRNA) %in% rownames(cds),]
  sumOverlap <- dim(subRNA)[1]
  if(sumOverlap == 0){
    stop(paste0("No overlaps between bulk RNA data and reduce dimensions feature found.",
                    "\nEither recreate counts matrix or most likely these data sets are incompatible!"))
  }
  if( (sumOverlap / rownames(cds)) < 0.25 ){
    if(force){
      warning("Less than 25% of the features are present in this bulk RNA data set! Continuing since force = TRUE!")
    }else{
      stop("Less than 25% of the features are present in this bulk RNA data set! Set force = TRUE to continue!")
    }
  }
  .logMessage("Overlap Ratio of Reduced Dims Features = ", (sumOverlap / length(rDGR)), verbose = TRUE, logFile = logFile)
  
  ##################################################
  # Create Bulk Matrix
  ##################################################
  bulkMat <- safeSubset(
    mat = getAssay(subRNA), 
    subsetRows = rownames(cds))

  ##################################################
  # Simulate and Project
  ##################################################
  depthN <- round(sum(rowSums(counts(cds)) / dim(cds)[2]))
  nRep <- 5
  n2 <- ceiling(n / nRep)
  ratios <- c(2, 1.5, 1, 0.5, 0.25) #range of ratios of number of fragments
  x<-1
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
    simRD <- as.matrix(.projectLSI(simMat, LSI = rD, verbose = verbose, binarize = binarize, LSIMethod=LSIMethod))
    rownames(simRD) <- paste0(colnames(bulkMat)[x], "#", seq_len(nrow(simRD)))
    simRD
  }, mc.cores =  threads) %>% Reduce("rbind", .)
  
  if(is.null(embedding)){
    if(rD$scaleDims){
      simRD <- .scaleDims(simRD)
    }
    out <- SimpleList(
      simulatedReducedDims = simRD
    )
    return(out)
  }
  .logThis(simRD, "simulatedReducedDims", logFile = logFile)
  
  ##################################################
  # Prep Reduced Dims
  ##################################################
  embedding <- getEmbedding(cds = cds, embedding = embedding, returnDF = FALSE)
  corCutOff <- embedding$params$corCutOff
  dimsToUse <- embedding$params$dimsToUse
  scaleDims <- embedding$params$scaleDims
  
  if(is.null(scaleDims)){
    scaleDims <- rD$scaleDims
  }
  
  simRD <- .scaleDims(simRD)
  
  if(embedding$params$nc != ncol(simRD)){
    
    if(is.null(dimsToUse)){
      dimsToUse <- seq_len(ncol(rD[[1]]))
    }
    
    if(!is.null(corCutOff)){
      if(scaleDims){
        corToDepth <- rD$corToDepth$scaled
        dimsToUse <- dimsToUse[corToDepth < corCutOff]
      }else{
        corToDepth <- rD$corToDepth$none
        dimsToUse <- dimsToUse[corToDepth < corCutOff]
      }
    }
    
    if(embedding$params$nc != ncol(simRD)){
      .logMessage("Error incosistency found with matching LSI dimensions to those used in addEmbedding",
                  "\nReturning with simulated reduced dimension coordinates...", verbose = TRUE, logFile = logFile)
      out <- SimpleList(
        simulatedReducedDims = simRD
      )
      return(out)
    }
    
    simRD <- simRD[, dimsToUse, drop = FALSE]
    
  }
  
  ##################################################
  # Get Previous UMAP Model
  ##################################################
  umapModel <- .loadUWOT(embedding$params$uwotModel, embedding$params$nc)
  
  idx <- sort(sample(seq_len(nrow(rD[[1]])), min(nrow(rD[[1]]), 5000))) #Try to use 5000 or total cells to check validity
  rD2 <- getReducedDims(
    cds = cds, 
    reducedDims = reducedDims, 
    dimsToUse = embedding$params$dimsToUse,
    scaleDims = embedding$params$scaleDims,
    corCutOff = embedding$params$corCutOff
  )[idx,,drop=FALSE]
  
  ##################################################
  # Project UMAP
  ##################################################
  set.seed(1)
  threads2 <- max(floor(threads/2), 1)
  simUMAP <- uwot::umap_transform(
    X = rbind(rD2, simRD), 
    model = umapModel, 
    verbose = TRUE, 
    n_threads = threads2
  )
  rownames(simUMAP) <- c(rownames(rD2), rownames(simRD))
  .logThis(simUMAP, "simulatedUMAP", logFile = logFile)
  
  #Check if the projection matches using previous data
  c1 <- cor(simUMAP[rownames(rD2), 1], embedding[[1]][rownames(rD2),1])
  c2 <- cor(simUMAP[rownames(rD2), 2], embedding[[1]][rownames(rD2),2])
  if(min(c1, c2) < 0.8){
    .logMessage("Warning projection correlation is less than 0.8 (R = ", round(min(c1,c2), 4),").\nThese results may not be accurate because of the lack of heterogeneity in the single cell data.", verbose = TRUE, logFile = logFile)
  }
  
  dfUMAP <- embedding[[1]]
  colnames(dfUMAP) <- c("UMAP1", "UMAP2")
  colnames(simUMAP) <- c("UMAP1", "UMAP2")
  dfUMAP <- DataFrame(dfUMAP)
  dfUMAP$Type <- Rle("scATAC", lengths = nrow(dfUMAP))
  
  simUMAP <- DataFrame(simUMAP[rownames(simRD),,drop=FALSE])
  simUMAP$Type <- Rle(stringr::str_split(rownames(simUMAP), pattern = "#", simplify = TRUE)[,1])
  
  out <- SimpleList(
    simulatedBulkUMAP = simUMAP,
    singleCellUMAP = dfUMAP,
    simulatedReducedDims = simRD
  )
  .endLogging(logFile = logFile)
  
  return(out)
  
}

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

safeSubset<-function (mat = NULL, subsetRows = NULL, subsetCols = NULL) 
{
  if (!is.null(subsetRows)) {
    idxNotIn <- which(!subsetRows %in% rownames(mat))
    if (length(idxNotIn) > 0) {
      subsetNamesNotIn <- subsetRows[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, 
                                       dims = c(length(idxNotIn), ncol = ncol(mat)))
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

projectLSI<-function (mat = NULL, LSI = NULL, LSIMethod=2, binarize=FALSE, returnModel = FALSE, verbose = FALSE) 
{
  .logThis(append(args, mget(names(formals()), sys.frame(sys.nframe()))), 
           "LSI-Projection Parameters", logFile = logFile)
  out2 <- tryCatch({
    require(Matrix)
    set.seed(LSI$seed)
    if (is.null(tstart)) {
      tstart <- Sys.time()
    }
    if(verbose) message(sprintf("Projecting LSI, Input Matrix = %s GB", 
                         round(object.size(mat)/10^9, 3)))
    if(verbose) message("Subsetting by Non-Zero features in inital Matrix")
    #mat <- mat[LSI$idx, ] I don't see what this does after safeSubset
    if (binarize) {
      .logDiffTime("Binarizing Matrix", tstart, addHeader = FALSE, 
                   verbose = verbose, logFile = logFile)
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
    if (LSIMethod == 1) {
      .if(verbose) message("Computing Inverse Document Frequency")
      idf <- as(log(1 + LSI$nCol/LSI$rowSm), "sparseVector")
      .if(verbose) message("Computing TF-IDF Matrix")
      mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*% 
        mat
    }
    else if (LSIMethod == 2) {
      if(verbose) message("Computing Inverse Document Frequency")
      idf <- as(LSI$nCol/LSI$rowSm, "sparseVector")
      .if(verbose) message("Computing TF-IDF Matrix")
      mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*% 
        mat
      mat@x <- log(mat@x * LSI$scaleTo + 1)
    }
    else if (LSIMethod == 3) {
      mat@x <- log(mat@x + 1)
      if(verbose) message("Computing Inverse Document Frequency")
      idf <- as(log(1 + LSI$nCol/LSI$rowSm), "sparseVector")
      if(verbose) message("Computing TF-IDF Matrix")
      mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*% 
        mat
    }
    else {
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
    .logDiffTime("Computing Projected Coordinates", tstart, 
                 addHeader = FALSE, verbose = verbose, logFile = logFile)
    svdDiag <- matrix(0, nrow = LSI$nDimensions, ncol = LSI$nDimensions)
    diag(svdDiag) <- LSI$svd$d
    matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
    matSVD <- as.matrix(matSVD)
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI", seq_len(ncol(matSVD)))
    if (returnModel) {
      .logDiffTime("Calculating Re-Projected Matrix", tstart, 
                   addHeader = FALSE, verbose = verbose, logFile = logFile)
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
    .logError(e, fn = ".projectLSI", info = "", errorList = errorList, 
              logFile = logFile)
  })
  out2
}
