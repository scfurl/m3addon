writeH5AD<-function (object, file, assay = NULL, graph = NULL, verbose = TRUE, 
          overwrite = FALSE, x.slot="counts", ...) 
{
  #message("WriteH5AD is not currently operational, please use as.loom")
  #.NotYetImplemented()
  # if (!PackageCheck("hdf5r", error = FALSE)) {
  #   stop("Please install hdf5r to enable h5ad functionality")
  # }
  # Seurat:::CheckDots(...)
  # if (file.exists(file) && !overwrite) {
  #   stop("Output file exists, not overwriting")
  # }
  # assay <- assay %||% DefaultAssay(object = object)
  # graph <- graph %||% paste0(assay, "_snn")
  # DefaultAssay(object = object) <- assay
  # object[["active_assay"]] <- Idents(object = object)
  # x.slot <- if (!IsMatrixEmpty(x = GetAssayData(object = object, 
  #                                               slot = "scale.data"))) {
  #   "scale.data"
  # }
  # else if (identical(x = GetAssayData(object = object, slot = "counts"), 
  #                    y = GetAssayData(object = object, slot = "data"))) {
  #   "counts"
  # }
  # else {
  #   "data"
  # }
  # if (verbose) {
  #   message("Storing '", x.slot, "' into 'X'")
  # }
  # raw.slot <- switch(EXPR = x.slot, scale.data = "data", data = "counts", 
  #                    NULL)
  # if (verbose) {
  #   message("Storing '", raw.slot, "' into 'raw.X'")
  # }
  # if (x.slot == "counts" && IsMatrixEmpty(x = GetAssayData(object = object, 
  #                                                          slot = x.slot))) {
  #   if (verbose) {
  #     warning("Counts slot empty, storing data slot into 'X', not storing a raw.X")
  #   }
  #   x.slot <- "data"
  #   raw.slot <- NULL
  # }
  meta.features <- object[[assay]][[]]
  colnames(x = meta.features) <- gsub(pattern = "dispersion.scaled", 
                                      replacement = "dispersions_norm", x = colnames(x = meta.features))
  colnames(x = meta.features) <- gsub(pattern = "dispersion", 
                                      replacement = "dispersions", x = colnames(x = meta.features))
  colnames(x = meta.features) <- gsub(pattern = "mean", replacement = "means", 
                                      x = colnames(x = meta.features))
  colnames(x = meta.features) <- gsub(pattern = "\\.", replacement = "_", 
                                      x = colnames(x = meta.features))
  meta.features$highly_variable <- FALSE
  meta.features[VariableFeatures(object = object), "highly_variable"] <- TRUE
  meta.features$index <- rownames(x = meta.features)
  mf.order <- c("index", grep(pattern = "index", x = colnames(x = meta.features), 
                              invert = TRUE, value = TRUE))
  meta.features <- meta.features[, mf.order, drop = FALSE]
  meta.data <- object[[]]
  assays.remove <- grep(pattern = assay, x = Seurat:::FilterObjects(object = object, 
                                                           classes.keep = "Assay"), invert = TRUE, value = TRUE)
  if (length(x = assays.remove)) {
    assays.remove <- grep(pattern = assays.remove, x = colnames(x = meta.data))
    meta.data <- meta.data[, -assays.remove, drop = FALSE]
  }
  colnames(x = meta.data) <- gsub(pattern = paste0("nCount_", 
                                                   assay), replacement = "n_counts", x = colnames(x = meta.data))
  colnames(x = meta.data) <- gsub(pattern = paste0("nFeatures_", 
                                                   assay), replacement = "n_umis", x = colnames(x = meta.data))
  colnames(x = meta.data) <- gsub(pattern = "\\.", replacement = "_", 
                                  x = colnames(x = meta.data))
  meta.data$index <- rownames(x = meta.data)
  md.order <- c("index", grep(pattern = "index", x = colnames(x = meta.data), 
                              invert = TRUE, value = TRUE))
  meta.data <- meta.data[, md.order, drop = FALSE]
  hfile <- hdf5r::h5file(filename = file, mode = "w")
  if (verbose) {
    message("Writing 'X' matrix")
  }
  x.data <- GetAssayData(object = object, slot = x.slot, assay = assay)
  switch(EXPR = x.slot, scale.data = hfile[["X"]] <- x.data, 
         {
           x.data <- as.sparse(x = x.data)
           hfile[["X/indices"]] <- slot(object = x.data, "i") - 
             1
           hfile[["X/indptr"]] <- slot(object = x.data, "p")
           hfile[["X/data"]] <- slot(object = x.data, "x")
         })
  if (verbose) {
    message("Writing 'var' metadata")
  }
  hfile[["var"]] <- meta.features[rownames(x = x.data), , drop = FALSE]
  if (!is.null(x = raw.slot)) {
    if (verbose) {
      message("Writing 'raw.X' sparse matrix")
    }
    raw.data <- GetAssayData(object = object, slot = raw.slot, 
                             assay = assay)
    hfile[["raw.X/indices"]] <- slot(object = raw.data, "i") - 
      1
    hfile[["raw.X/indptr"]] <- slot(object = raw.data, "p")
    hfile[["raw.X/data"]] <- slot(object = raw.data, "x")
    if (verbose) {
      message("Writing 'raw.var' metadata")
    }
    hfile[["raw.var"]] <- meta.features
  }
  if (verbose) {
    message("Writing 'obs' metadata")
  }
  hfile[["obs"]] <- meta.data
  if (x.slot == "scale.data") {
    dim.reducs <- FilterObjects(object = object, classes.keep = "DimReduc")
    dim.reducs <- Filter(f = function(x) {
      return(DefaultAssay(object = object[[x]]) == assay)
    }, x = dim.reducs)
    if (length(x = dim.reducs) >= 1) {
      embedding.names <- paste0("X_", dim.reducs)
      names(x = embedding.names) <- dim.reducs
      loading.names <- gsub(pattern = "_$", replacement = "s", 
                            x = vapply(X = dim.reducs, FUN = function(x) {
                              return(Key(object[[x]]))
                            }, FUN.VALUE = character(length = 1L)))
      embeddings <- sapply(X = dim.reducs, FUN = function(x) {
        return(t(x = Embeddings(object = object, reduction = x)))
      }, USE.NAMES = TRUE, simplify = FALSE)
      names(x = embeddings) <- embedding.names[names(x = embeddings)]
      hfile[["obsm"]] <- embeddings
    }
    else if (verbose) {
      warning("No dimensional reduction objects for assay '", 
              assay, "' found")
    }
  }
  else if (verbose) {
    warning("Intial object unscaled, not storing dimensional reduction information")
  }
  if (x.slot == "scale.data") {
    graphs <- FilterObjects(object = object, classes.keep = "Graph")
    graphs <- grep(pattern = graph, x = graphs, value = TRUE)
    if (length(x = graphs) == 1) {
      ""
    }
    else if (verbose) {
      warning("Could not find a graph named '", graph, 
              "'")
    }
  }
  else if (verbose) {
    warning("Initial object unscaled, not storing graph information")
  }
  hfile$flush()
  invisible(x = hfile)
}
