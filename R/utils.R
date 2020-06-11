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
#' @param cds_list Input cell_data_set object or sparse matrix.
#' @import Matrix
#' @export
tf_idf_transform <- function(input){
  if(class(input)=="cell_data_set"){
    mat<-exprs(input)
  }else{
    mat<-input
  }
  idf <- log( ncol(mat) / ( 1 + Matrix::rowSums(mat != 0) ) )
  idf<-.sparseDiagonal(x=idf)
  tf_idf <- crossprod(mat, idf)
  colnames(tf_idf) <- rownames(mat)
  tf_idf_out<-Matrix::t(tf_idf / sqrt( Matrix::rowSums( tf_idf^2 ) ))
  if(class(input)=="cell_data_set"){
    input@assays$data$counts<-tf_idf_out
    return(input)
  }else{
    return(tf_idf_out)
  }
}

#' Performs TF-IDF transformation on a cell_data_set v2
#'
#' @description Just like it sounds but different.
#'
#' @param cds_list Input cell_data_set object or sparse matrix.
#' @import Matrix
#' @export
tf_idf_transform_v2 <- function(input){
  if(class(input)=="cell_data_set"){
    mat<-exprs(input)
  }else{
    mat<-input
  }
  colSm <- Matrix::colSums(mat)
  rowSm <- Matrix::rowSums(mat)
  freqs <- t(t(mat)/colSm)
  idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  tfidf@x[is.na(tfidf@x)] <- 0
  if(class(input)=="cell_data_set"){
    input@assays$data$counts<-tfidf
    return(input)
  }else{
    return(tfidf)
  }
}

#' @export
svd_lsi<-function(sp_mat, num_dim, mat_only=T){
  svd <- irlba::irlba(sp_mat, num_dim, num_dim)
  svdDiag <- matrix(0, nrow=num_dim, ncol=num_dim)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(sp_mat)
  colnames(matSVD) <- seq_len(ncol(matSVD))
  if(mat_only){
    return(matSVD)
  }else{
    return(list(matSVD=matSVD, svd=svd))
  }
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

#' @export
bimodality_coefficient<-function(x, finite=TRUE,...){
  if(finite==TRUE){
    G=skewness(x,finite)
    sample.excess.kurtosis=kurtosis(x,finite)
    K=sample.excess.kurtosis
    n=length(x)
    B=((G^2)+1)/(K+ ((3*((n-1)^2))/((n-2)*(n-3))))
  }
  else{
    G=skewness(x,FALSE)
    K=kurtosis(x,FALSE)
    B=((G^2)+1)/(K)
  }
  return(B)
}

#' @export
skewness<-function(x, finite=TRUE){
  n=length(x)
  S=(1/n)*sum((x-mean(x))^3)/(((1/n)*sum((x-mean(x))^2))^1.5)
  if(finite==FALSE){
    S=S
  }else{
    S=S*(sqrt(n*(n-1)))/(n-2)
  }
  return(S)	
}

#' @export
kurtosis<-function(x, finite){
  n=length(x)
  K=(1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2) - 3
  if(finite==FALSE){
    K=K
  }
  else{
    K=((n-1)*((n+1)*K - 3*(n-1))/((n-2)*(n-3))) +3
  }
  return(K)	
}

#' @export

lighten_darken_color<-function(col, amt) {
  if (substring(col, 1, 1)=="#") {
    col = substring(col, 2)
  }
  num = as.hexmode(col)
  r = bitwShiftR(num, 16) + amt
  if (r > 255) {r = 255}
  if  (r < 0) {r = 0}
  b = bitwAnd(bitwShiftR(num, 8), 0x00FF) + amt
  if (b > 255) {b = 255}
  if  (b < 0) {b = 0}
  g = bitwAnd(num, 0x0000FF) + amt
  if (g > 255) {g = 255}
  if (g < 0) {g = 0}
  inter<-paste("000000", as.hexmode(bitwOr(g , bitwOr(bitwShiftL(b, 8), bitwShiftL(r, 16)))), sep="")
  ret<-substr(inter, nchar(inter)-5, nchar(inter))
  return(toupper(paste("#", ret, sep="")))
}




