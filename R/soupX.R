#' @title Load Data for SoupX
#' @description Loads unfiltered 10X data from each data-set and identifies which droplets are cells using the cellranger defaults.
#' @param dataDir Top level cellranger output directory (the directory that contains the "raw_gene_bc_matrices" folder).
#' @param cellIDs	Barcodes of droplets that contain cells. If NULL, use the default cellranger set.
#' @param channelName The name of the channel to store. If NULL set to either names(dataDir) or dataDir is no name is set.
#' @param readArgs A list of extra parameters passed to Seurat::Read10X.
#' @param includeFeatures If multiple feature types are present, keep only the types mentioned here and collapse to a single matrix.
#' @param method type of input data; Use "SoupX_10X" for processing cellranger output, "10X_H5" for processing a folder of H5 
#' cellranger output files, and "STARsolo" for processing STARsolo output.
#' @param verbose Be verbose?
#' @param ... Extra parameters passed to SoupChannel construction function.
#' @return A SoupChannel object containing the count tables for the 10X dataset.
#' @importFrom Seurat Read10X
#' @export
#dataDir<-"/Users/sfurla/OneDrive - Fred Hutchinson Cancer Research Center/computation/Analysis/BMME/KSTd145/data/STARsolo/L1_SO/outs/Gene"
loadSoupX<-function (dataDir, cellIDs = NULL, channelName = NULL, readArgs = list(),
          includeFeatures = c("Gene Expression"), verbose = TRUE, method=c("SoupX_10X", "10X_H5", "STARsolo"), ...)
{
  if (verbose) 
    message(sprintf("Loading raw count data"))
  if(method=="SoupX_10X"){
    isV3 = dir.exists(file.path(dataDir, "raw_feature_bc_matrix"))
    tgt = file.path(dataDir, ifelse(isV3, "raw_feature_bc_matrix", 
                                    "raw_gene_bc_matrices"))
    if (!isV3) 
      tgt = file.path(tgt, list.files(tgt))
    dat = do.call(Read10X, c(list(data.dir = tgt), readArgs))
  }
  if(method=="STARsolo"){
    tgt = file.path(dataDir, "raw")
    dat = do.call(readSTARsolo, c(list(data.dir = tgt), readArgs))
  }
  if(method=="10X_H5"){
    isV3 = file.exists(file.path(dataDir, "raw_feature_bc_matrix.h5"))
    tgt = file.path(dataDir, ifelse(isV3, "raw_feature_bc_matrix.h5", 
                                    "raw_gene_bc_matrices.h5"))
    if (!isV3) 
      tgt = file.path(tgt, list.files(tgt))
    dat = read.cellranger.h5.file(tgt)
  }
  if (verbose) 
    message(sprintf("Loading cell-only count data"))
  if (!is.null(cellIDs)) {
    if (all(grepl("\\-1$", cellIDs))) 
      cellIDs = gsub("\\-1$", "", cellIDs)
    if (!all(cellIDs %in% colnames(dat))) 
      stop("Not all supplied cellIDs found in raw data.")
    datCells = dat[, match(cellIDs, colnames(dat))]
  }else {
    if(method=="SoupX_10X"){
      tgt = file.path(dataDir, ifelse(isV3, "filtered_feature_bc_matrix", 
                                      "filtered_gene_bc_matrices"))
      if (!isV3) 
        tgt = file.path(tgt, list.files(tgt))
      datCells = do.call(Read10X, c(list(data.dir = tgt), readArgs))
      if (is.list(dat)) {
        dat = do.call(rbind, dat[includeFeatures])
        datCells = do.call(rbind, datCells[includeFeatures])
      }
    }
    if(method=="10X_H5"){
      tgt = file.path(dataDir, ifelse(isV3, "filtered_feature_bc_matrix.h5", 
                                      "filtered_gene_bc_matrices.h5"))
      if (!isV3) {tgt = file.path(tgt, list.files(tgt))}
      datCells = read.cellranger.h5.file(tgt)
    }
    if(method=="STARsolo"){
      tgt = file.path(dataDir, "filtered")
      datCells = do.call(readSTARsolo, c(list(data.dir = tgt), readArgs))
    }
  }
  if (verbose) 
    message(sprintf("Loading extra analysis data where available"))
  mDat = NULL
  tgt = file.path(dataDir, "analysis", "clustering", "graphclust", 
                  "clusters.csv")
  if (file.exists(tgt)) {
    clusters = read.csv(tgt)
    mDat = data.frame(clusters = clusters$Cluster, row.names = clusters$Barcode)
  }
  tgt = file.path(dataDir, "analysis", "clustering", "kmeans_10_clusters", 
                  "clusters.csv")
  if (file.exists(tgt)) {
    clusters = read.csv(tgt)
    mDat$clustersFine = clusters$Cluster
  }
  tgt = file.path(dataDir, "analysis", "tsne", "2_components", 
                  "projection.csv")
  if (file.exists(tgt)) {
    tsne = read.csv(tgt)
    if (is.null(mDat)) {
      mDat = data.frame(tSNE1 = tsne$TSNE.1, tSNE2 = tsne$TSNE.2, 
                        row.names = tsne$Barcode)
    }
    else {
      mDat$tSNE1 = tsne$TSNE.1[match(rownames(mDat), tsne$Barcode)]
      mDat$tSNE2 = tsne$TSNE.2[match(rownames(mDat), tsne$Barcode)]
    }
    DR = c("tSNE1", "tSNE2")
  }
  else {
    DR = NULL
  }
  if (!is.null(mDat) && any(rownames(mDat) != colnames(datCells))) {
    rownames(mDat) = gsub("-1$", "", rownames(mDat))
    if (any(rownames(mDat) != colnames(datCells))) 
      stop("Error matching meta-data to cell names.")
  }
  if (is.null(channelName)) 
    channelName = ifelse(is.null(names(dataDir)), dataDir, 
                         names(dataDir))
  channel = SoupChannel(tod = dat, toc = datCells, metaData = mDat, 
                        channelName = channelName, dataDir = dataDir, dataType = "10X", 
                        isV3 = isV3, DR = DR, ...)
  return(channel)
}




#' @import hdf5r
#' @importFrom Matrix sparseMatrix
read.cellranger.h5.file = function(h5.file, gene.column = "gene_names", unique.features = TRUE) {
  if(!file.exists(h5.file)){stop(paste0("File ", h5.file," not found"))}
  #s<-h5file(h5.file, mode="r")
  s <- H5File$new(h5.file, mode="r+")
  
  #s$ls(recursive=T)
  genome=s[["matrix/features/genome"]][]
  barcodes = s[["matrix/barcodes"]][]
  gene_ids = s[["matrix/features/id"]][]
  gene_names =s[["matrix/features/name"]][]
  featuretype = s[["matrix/features/feature_type"]][]
  data = s[["matrix/data"]][]
  indices = s[["matrix/indices"]][]+1
  indptr = s[["matrix/indptr"]][]
  shape = s[["matrix/shape"]][]
  #browser()
  h5close(s)
  # gbm = new(
  #   "dgCMatrix",
  #   x = data, i = indices, p = indptr,
  #   Dim = shape)
  if(length(levels(factor(featuretype)))>1){
    warning("Warning - Multiple feature types found")
  }
  
  gbm = sparseMatrix(x = data, i = indices, p = indptr,
                     dims = shape)
  
  colnames(gbm) = barcodes
  # pre_ver_3<-F
  # features.loc<-file.path(tgt, "features.tsv.gz")
  # feature.names1 <- read.delim(file = ifelse(test = pre_ver_3, 
  #                                           yes = gene.loc, no = features.loc), header = FALSE, 
  #                             stringsAsFactors = FALSE)
  feature.names<-data.frame(gene_ids, gene_names, featuretype, stringsAsFactors = F)
  if (any(is.na(x = feature.names[, gene.column]))) {
    warning("Some features names are NA. Replacing NA names with ID from the opposite column requested", 
            call. = FALSE, immediate. = TRUE)
    na.features <- which(x = is.na(x = feature.names[, 
                                                     gene.column]))
    replacement.column <- ifelse(test = gene.column == 
                                   2, yes = 1, no = 2)
    feature.names[na.features, gene.column] <- feature.names[na.features, 
                                                             replacement.column]
  }
  if (unique.features) {
    fcols = ncol(x = feature.names)
    if (fcols < match(gene.column, colnames(feature.names))) {
      stop(paste0("gene.column was set to ", gene.column, 
                  " but feature.tsv.gz (or genes.tsv) only has ", 
                  fcols, " columns.", " Try setting the gene.column argument to a value <= to ", 
                  fcols, "."))
    }
    rownames(x = gbm) <- make.unique(names = feature.names[, 
                                                            gene.column])
  }
  if (ncol(x = feature.names) > 2) {
    data_types <- factor(x = feature.names$featuretype)
    lvls <- levels(x = data_types)
    if (length(x = lvls) > 1) {
      message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
    }
    expr_name <- "Gene Expression"
    if (expr_name %in% lvls) {
      lvls <- c(expr_name, lvls[-which(x = lvls == 
                                         expr_name)])
    }
    gbm <- lapply(X = lvls, FUN = function(l) {
      return(gbm[data_types == l, ])
    })
    names(x = gbm) <- lvls
  }
  else {
    gbm <- list(gbm)
  }
  gbm
}


readSTARsolo<-function (data.dir = NULL, gene.column = 2, unique.features = TRUE) 
{
  full.data <- list()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, "barcodes.tsv")
    features.loc <- file.path(run, "features.tsv")
    matrix.loc <- file.path(run, "matrix.mtx")
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc)) {
      stop("Gene name or features file missing. Expecting ", 
           basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", 
           basename(path = matrix.loc))
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    if (is.null(x = names(x = data.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], 
                                   "_", cell.names)
    }
    feature.names <- read.delim(file = ifelse(test = pre_ver_3, 
                                              yes = gene.loc, no = features.loc), header = FALSE, 
                                stringsAsFactors = FALSE)
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning("Some features names are NA. Replacing NA names with ID from the opposite column requested", 
              call. = FALSE, immediate. = TRUE)
      na.features <- which(x = is.na(x = feature.names[, 
                                                       gene.column]))
      replacement.column <- ifelse(test = gene.column == 
                                     2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, 
                                                               replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column, 
                    " but feature.tsv.gz (or genes.tsv) only has ", 
                    fcols, " columns.", " Try setting the gene.column argument to a value <= to ", 
                    fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, 
                                                              gene.column])
    }
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 
          0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls == 
                                           expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, ])
      })
      names(x = data) <- lvls
    }
    else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, 
                                               FUN = `[[`, j))
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  }
  else {
    return(list_of_data)
  }
}
