#' Finds common features in a list of cds objects
#'
#' @description Machine learning algorithms often require features to be the same across 
#' datasets.  Thisfunction finds common features between a list of cell data set objects and 
#' returns a list of cds's that have the same features.  Note that this function uses rownames 
#' of the 'fData' DataFrame to find the intersect of features common to all cds's
#'
#' @param cds_list Input cell_data_set object.
#' @export
common_features <- function(cds_list){
  len<-length(cds_list)
  common_features=vector()
  for(i in 1:len){
    if(i < 2){
      common_features<-rownames(fData(cds_list[[i]]))
    }else{
      common_features<-unique(intersect(common_features, rownames(fData(cds_list[[i]]))))
    }
  }
  for(i in 1:len){
    cds_list[[i]]<-cds_list[[i]][match(common_features, rownames(cds_list[[i]])),]
  }
  return(cds_list)
}



#' Performs TF-IDF transformation on a cell_data_set
#'
#' @description Just like it sounds.
#'
#' @param cds_list Input cell_data_set object.
#' @import Matrix
#' @export
tf_idf_transform <- function(cds){
  idf <- log( ncol(exprs(cds)) / ( 1 + Matrix::rowSums(exprs(cds) != 0) ) )
  idf<-.sparseDiagonal(x=idf)
  tf_idf <- crossprod(exprs(cds), idf)
  colnames(tf_idf) <- rownames(exprs(cds))
  tf_idf_out<-Matrix::t(tf_idf / sqrt( Matrix::rowSums( tf_idf^2 ) ))
  cds@assays$data$counts<-tf_idf_out
  return(cds)
}


Noisify <- function(data, amount=0.0001) {
  if (is.vector(data)) {
    noise <- runif(length(data), -amount, amount)
    noisified <- data + noise
  } else {
    length <- dim(data)[1] * dim(data)[2]
    noise <- matrix(runif(length, -amount, amount), dim(data)[1])
    noisified <- data + noise
  }
  return(noisified)
}


#' Detects genes above minimum threshold.
#'
#' @description For each gene in a cell_data_set object, detect_genes counts
#' how many cells are expressed above a minimum threshold. In addition, for
#' each cell, detect_genes counts the number of genes above this threshold that
#' are detectable. Results are added as columns num_cells_expressed and
#' num_genes_expressed in the rowData and colData tables respectively.
#'
#' @param cds Input cell_data_set object.
#' @param min_expr Numeric indicating expression threshold
#' @param exprs_bin Boolean whether to bin genes by mean expression
#' @param exprs_cuts Numeic indicating number of bins if using exprs_bin
#' @return Updated cell_data_set object
#' @importFrom Hmisc cut2
#' @export
detect_genes <- function(cds, min_expr=0, exprs_bin=TRUE, exprs_cuts=25){
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(is.numeric(min_expr))
  
  rowData(cds)$num_cells_expressed <- Matrix::rowSums(SingleCellExperiment::counts(cds) > min_expr)
  colData(cds)$num_genes_expressed <- Matrix::colSums(SingleCellExperiment::counts(cds) > min_expr)
  if(exprs_bin){
    fData(cds)$exprs_bin = cut2(log(Matrix::rowMeans(normalized_counts(cds))), m=floor(nrow(fData(cds))/exprs_cuts))
  }
  cds
}


#' Convert a GMT File to a ragged list of genes
#'
#' @description GMT Files (See MSigDB) are a convenient means of storing genesets
#'
#' @param GMTfn GMT filename
#' @export
GMT_to_list<-function (GMTfn) 
{
  file.data <- readLines(GMTfn)
  num.lines <- length(file.data)
  output <- vector("list", num.lines)
  for (i in 1:num.lines) {
    vec <- unlist(strsplit(file.data[i], "\t"))
    output[[i]] <- vec[3:length(vec)]
    names(output)[i] <- vec[1]
  }
  return(output)
}
