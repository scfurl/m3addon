#' 
#' #' @param standardize Standardize matrices - scales columns to have unit variance
#' #' and mean 0
#' #' @param num.cc Number of canonical vectors to calculate
#' #' @param seed.use Random seed to set. If NULL, does not set a seed
#' #' @param verbose Show progress messages
#' #'
#' #' @importFrom irlba irlba
#' #'
#' #' @rdname RunCCA
#' #' @concept dimensional_reduction
#' #' @export
#' #'
#' RunCCA_engine <- function(
#'   object1,
#'   object2,
#'   standardize = F,
#'   num.cc = 20,
#'   seed.use = 42,
#'   verbose = FALSE,
#'   ...
#' ) {
#'   if (!is.null(x = seed.use)) {
#'     set.seed(seed = seed.use)
#'   }
#'   cells1 <- colnames(x = object1)
#'   cells2 <- colnames(x = object2)
#'   ##fix later
#'   if (standardize) {
#'     object1 <- Standardize(mat = object1, display_progress = FALSE)
#'     object2 <- Standardize(mat = object2, display_progress = FALSE)
#'   }
#'   mat3 <- crossprod(x = object1, y = object2)
#'   cca.svd <- irlba(A = mat3, nv = num.cc)
#'   cca.data <- rbind(cca.svd$u, cca.svd$v)
#'   colnames(x = cca.data) <- paste0("CC", 1:num.cc)
#'   rownames(cca.data) <- c(cells1, cells2)
#'   cca.data <- apply(
#'     X = cca.data,
#'     MARGIN = 2,
#'     FUN = function(x) {
#'       if (sign(x[1]) == -1) {
#'         x <- x * -1
#'       }
#'       return(x)
#'     }
#'   )
#'   return(list(ccv = cca.data, d = cca.svd$d))
#' }
#' 
#' #' @title RunCCA
#' #' @param assay1,assay2 Assays to pull from in the first and second objects, respectively
#' #' @param features Set of genes to use in CCA. Default is the union of both
#' #' the variable features sets present in both objects.
#' #' @param renormalize Renormalize raw data after merging the objects. If FALSE,
#' #' merge the data matrices also.
#' #' @param rescale Rescale the datasets prior to CCA. If FALSE, uses existing data in the scale data slots.
#' #' @param compute.gene.loadings Also compute the gene loadings. NOTE - this will
#' #' scale every gene in the dataset which may impose a high memory cost.
#' #' @param add.cell.id1,add.cell.id2 Add ...
#' #' @param ... Extra parameters (passed onto MergeSeurat in case with two objects
#' #' passed, passed onto ScaleData in case with single object and rescale.groups
#' #' set to TRUE)
#' #' @rdname RunCCA
#' #' @concept dimensional_reduction
#' #' @export
#' #'
#' RunCCA <- function(
#'   object1,
#'   object2,
#'   assay1 = NULL,
#'   assay2 = NULL,
#'   num.cc = 20,
#'   features = NULL,
#'   renormalize = FALSE,
#'   rescale = FALSE,
#'   compute.gene.loadings = TRUE,
#'   add.cell.id1 = NULL,
#'   add.cell.id2 = NULL,
#'   verbose = TRUE,
#'   ...
#' ) {
#' 
#'   if (is.null(x = features)) {
#'     if (length(x = VariableFeatures(object = object1, assay = assay1)) == 0) {
#'       stop(paste0("VariableFeatures not computed for the ", assay1, " assay in object1"))
#'     }
#'     if (length(x = VariableFeatures(object = object2, assay = assay2)) == 0) {
#'       stop(paste0("VariableFeatures not computed for the ", assay2, " assay in object2"))
#'     }
#'     features <- union(x = VariableFeatures(object = object1), y = VariableFeatures(object = object2))
#'     if (length(x = features) == 0) {
#'       stop("Zero features in the union of the VariableFeature sets ")
#'     }
#'   }
#'   nfeatures <- length(x = features)
#'   if (!(rescale)) {
#'     data.use1 <- GetAssayData(object = object1, assay = assay1, slot = "scale.data")
#'     data.use2 <- GetAssayData(object = object2, assay = assay2, slot = "scale.data")
#'     features <- CheckFeatures(data.use = data.use1, features = features, object.name = "object1", verbose = FALSE)
#'     features <- CheckFeatures(data.use = data.use2, features = features, object.name = "object2", verbose = FALSE)
#'     data1 <- data.use1[features, ]
#'     data2 <- data.use2[features, ]
#'   }
#'   if (rescale) {
#'     data.use1 <- GetAssayData(object = object1, assay = assay1, slot = "data")
#'     data.use2 <- GetAssayData(object = object2, assay = assay2, slot = "data")
#'     features <- CheckFeatures(data.use = data.use1, features = features, object.name = "object1", verbose = FALSE)
#'     features <- CheckFeatures(data.use = data.use2, features = features, object.name = "object2", verbose = FALSE)
#'     data1 <- data.use1[features,]
#'     data2 <- data.use2[features,]
#'     if (verbose) message("Rescaling groups")
#'     data1 <- FastRowScale(as.matrix(data1))
#'     dimnames(data1) <- list(features, colnames(x = object1))
#'     data2 <- FastRowScale(as.matrix(data2))
#'     dimnames(data2) <- list(features, colnames(x = object2))
#'   }
#'   if (length(x = features) / nfeatures < 0.1 & verbose) {
#'     warning("More than 10% of provided features filtered out. Please check that the given features are present in the scale.data slot for both the assays provided here and that they have non-zero variance.")
#'   }
#'   if (length(x = features) < 50) {
#'     warning("Fewer than 50 features used as input for CCA.")
#'   }
#'   if (verbose) {
#'     message("Running CCA")
#'   }
#'   cca.results <- RunCCA_engine(
#'     object1 = data1,
#'     object2 = data2,
#'     standardize = TRUE,
#'     num.cc = num.cc,
#'     verbose = verbose,
#'   )
#'   if (verbose) {
#'     message("Merging objects")
#'   }
#'   combined.object <- merge(
#'     x = object1,
#'     y = object2,
#'     merge.data = TRUE,
#'     ...
#'   )
#'   combined.object[['cca']] <- CreateDimReducObject(
#'     embeddings = cca.results$ccv[colnames(combined.object), ],
#'     assay = assay1,
#'     key = "CC_"
#'   )
#'   combined.object[['cca']]@assay.used <- DefaultAssay(combined.object)
#'   if (ncol(combined.object) != (ncol(object1) + ncol(object2))) {
#'     warning("Some cells removed after object merge due to minimum feature count cutoff")
#'   }
#'   combined.scale <- cbind(data1,data2)
#'   combined.object <- SetAssayData(object = combined.object,new.data = combined.scale, slot = "scale.data")
#'   if (renormalize) {
#'     combined.object <- NormalizeData(
#'       object = combined.object,
#'       assay = assay1,
#'       normalization.method = object1[[paste0("NormalizeData.", assay1)]]$normalization.method,
#'       scale.factor = object1[[paste0("NormalizeData.", assay1)]]$scale.factor
#'     )
#'   }
#'   if (compute.gene.loadings) {
#'     combined.object <- ProjectDim(
#'       object = combined.object,
#'       reduction = "cca",
#'       verbose = FALSE,
#'       overwrite = TRUE)
#'   }
#'   return(combined.object)
#' }
