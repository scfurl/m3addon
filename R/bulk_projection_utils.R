#' Extract data using a set of features
#' 
#' @description This function extracts data from a subject object (a SummarizedExperiment-like 
#' object) using the features specified in query.  The function is parameterized to allow for finding features by either feature
#' name (annotation-based) or using GenomicRanges (range-based).  Annotation-based searching defaults
#' to the rownames of the subject, but this can be altered using annotation argument
#' which will search for these colnames in the subject's feature metadata.  The default behavior
#' of this function is to return a SummarizedExperiment object with the full feature-set of the query and thus
#' will add zeros for those features in the query not found in the subject. Set the full-output argument to TRUE
#' to modify this behavior
#' 
#' @param query a vector of feature names or a GrangesObject
#' @param subject a SummarizedExperiment based object
#' @param annotation = A string specifying the target feature search for oin subject; Default is "row.names"
#' @param duplicate_hits parameter for dealing with multiple subject hits in a range-based search.  Options are to select those
#' multiple hits based on the following: "max.mean", "max.var", "max.disp", "min.mean", "min.var", "min.disp". Default is "max.disp"
#' @param fill_output = A boolean indicating whether the function should return a subject object filled with zero data for 
#' features in query that were not found in the subject
#' @param ignore_strand Whether to ignore strand for range-based searches
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can 
#' be set to FALSE for a cleaner output.
#' @return a list cotaining 1) matrix of the extracted data, 2) the ratio of features
#' found in the subject to those present in the query, 3) those features not_found
#' @export
#' 
extract_data<-function(query, 
                           subject,
                          annotation="row.names",
                          verbose=T,
                          fill_output=T,
                          duplicate_hits="max.disp",
                          ignore_strand=T){
  ##################################################
  # Check Inputs
  ##################################################
  check_input(input = query, name = "query", valid = c("character", "Granges"))
  check_input(input = subject, name = "subject", valid = c("SummarizedExperiment"))
  check_input(input = annotation, name = "annotation", valid = c("NULL", "character"))
  check_input(input = verbose, name = "verbose", valid = c("boolean"))
  check_input(input = fill_output, name = "fill_output", valid = c("boolean"))
  
  ##################################################
  # parse duplicate_hits behavior
  ################################################## 
  FUN<-get(strsplit(duplicate_hits, "\\.")[[1]][1])
  Var1<-strsplit(duplicate_hits, "\\.")[[1]][2]
  
  if(class(query)=="GRanges") {
      #rename features as coordinates are no longer meaningful
      fs<-names(query)<-paste("f", seq_along(query))
      se_sub <- subsetByOverlaps( subject, query, ignore.strand = ignore_strand, type="any")
      #Remove dups
      hits<-data.table::as.data.table(findOverlaps(query, se_sub, ignore.strand = TRUE))
      dat<-as.matrix(get_assay(se_sub))
      hits$var<-rowVars(dat)[hits$subjectHits]
      hits$mean<-rowMeans(dat)[hits$subjectHits]
      hits$disp<-sqrt(hits$var)/hits$mean*100
      best_hits <- hits %>% group_by(queryHits) %>% filter(get(Var1) == FUN(get(Var1)))
      if(nrow(best_hits) == 0){
        stop(paste0("No overlap between query and subject found."))
      }
      mat<-dat[best_hits$subjectHits,]
      rownames(mat)<-names(query)[best_hits$queryHits]
      if(fill_output)  mat<-safe_subset(mat, subsetRows = fs)
      overlap=nrow(best_hits)/length(fs)
      not_found<-query[!best_hits$queryHits %in% 1:length(names(fs)),]
  }else{
      if(annotation=="row.names") {
        fidx <- which(rownames(subject) %in% query)
        if(length(fidx) == 0){
          stop(paste0("No overlap between query and subject found."))
        }
        fs<-query
        mat<-as.matrix(get_assay(subject))[fidx,]
        if(is.null(rownames(mat))){
          rownames(mat)<-rownames(subject)[fidx]
        }
        if(fill_output)  mat<-safe_subset(mat, subsetRows = fs)
        overlap<-length(fidx)/length(query)
        not_found<-query[!fidx %in% 1:length(query)]
      }else{
        stop("not_implemented")
      }
  }
  return(list(mat=mat, overlap=overlap, notfound=not_found))
}


#'check_input helper
#'
#' @description
#' Adapted from: Jeffrey M. Granja, M. Ryan Corces, Sarah E. Pierce, S. Tansu Bagdatli, Hani Choudhry, Howard Y. Chang, William J. Greenleaf
#' doi: https://doi.org/10.1101/2020.04.28.066498
#' @export
check_input<-function (input = NULL, name = NULL, valid = NULL) 
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

#' safe_subset helper
#'
#' @description
#' Adapted from: Jeffrey M. Granja, M. Ryan Corces, Sarah E. Pierce, S. Tansu Bagdatli, Hani Choudhry, Howard Y. Chang, William J. Greenleaf
#' doi: https://doi.org/10.1101/2020.04.28.066498
#' 
#' @export
safe_subset<-function (mat = NULL, subsetRows = NULL, subsetCols = NULL) 
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

#' get_assay helper
#'
#' @description
#' Adapted from: Jeffrey M. Granja, M. Ryan Corces, Sarah E. Pierce, S. Tansu Bagdatli, Hani Choudhry, Howard Y. Chang, William J. Greenleaf
#' doi: https://doi.org/10.1101/2020.04.28.066498
#' 
#' @export
get_assay<-function (se = NULL, assayName = NULL) 
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
