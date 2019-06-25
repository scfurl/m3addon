#' @export
plot_gene_dispersion<-function(cds){
  prd<-cds@int_metadata$dispersion
  g<-ggplot2::ggplot(prd, ggplot2::aes(x = log_mean, y = fit)) 
  if("use_for_ordering" %in% colnames(cds@int_metadata$dispersion)){
    g <- g + ggplot2::geom_point(data=prd, ggplot2::aes(x=log_mean, y=log_dispersion, color=use_for_ordering), alpha=0.4)
  }else{
    g <- g + ggplot2::geom_point(data=prd, ggplot2::aes( x=log_mean, y=log_dispersion), color="grey", alpha=0.4)
  }
  g<-g+
    ggplot2::theme_bw() +
    ggplot2::geom_line() +
    ggplot2::geom_smooth(data=prd, ggplot2::aes(ymin = lci, ymax = uci), stat = "identity")
  g
}


#' Select genes in a cell_data_set for dimensionality reduction
#'
#' @description Monocle3 aims to learn how cells transition through a
#' biological program of gene expression changes in an experiment. Each cell
#' can be viewed as a point in a high-dimensional space, where each dimension
#' describes the expression of a different gene. Identifying the program of
#' gene expression changes is equivalent to learning a \emph{trajectory} that
#' the cells follow through this space. However, the more dimensions there are
#' in the analysis, the harder the trajectory is to learn. Fortunately, many
#' genes typically co-vary with one another, and so the dimensionality of the
#' data can be reduced with a wide variety of different algorithms. Monocle3
#' provides two different algorithms for dimensionality reduction via
#' \code{reduce_dimensions} (UMAP and tSNE). The function
#' \code{select_features} is an optional step in the trajectory building
#' process before \code{preprocess_cds}.  After calculating dispersion for
#' a cell_data_set using the \code{calculate_gene_dispersion} function, the 
#' \code{select_genes} function allows the user to identify a set of genes
#' that will be used in downstream dimensionality reduction methods.
#'
#'
#' @param cds the cell_data_set upon which to perform this operation.
#' @param fit_min the minimum multiple of the dispersion fit calculation; default = 1
#' @param fit_max the maximum multiple of the dispersion fit calculation; default = Inf
#' @return an updated cell_data_set object with selected features 
#' @export

select_genes<-function(cds, fit_min=1, fit_max=Inf, logcv_ll=NULL, logcv_ul=NULL, logmean_ul=NULL, logmean_ll=NULL){
  df<-cds@int_metadata$dispersion
  df$ratio<-df$log_dispersion/df$fit
  cds@int_metadata$dispersion$use_for_ordering <- df$ratio > fit_min & df$ratio < fit_max
  cds
}

#' @export
get_ordering_genes<-function(cds, gene_column="id"){
  cds@int_metadata$dispersion[[gene_column]][cds@int_metadata$dispersion$use_for_ordering]
}

#' Calculate dispersion genes in a cell_data_set object
#'
#' @description Monocle3 aims to learn how cells transition through a
#' biological program of gene expression changes in an experiment. Each cell
#' can be viewed as a point in a high-dimensional space, where each dimension
#' describes the expression of a different gene. Identifying the program of
#' gene expression changes is equivalent to learning a \emph{trajectory} that
#' the cells follow through this space. However, the more dimensions there are
#' in the analysis, the harder the trajectory is to learn. Fortunately, many
#' genes typically co-vary with one another, and so the dimensionality of the
#' data can be reduced with a wide variety of different algorithms. Monocle3
#' provides two different algorithms for dimensionality reduction via
#' \code{reduce_dimensions} (UMAP and tSNE). The function
#' \code{calculate_dispersion} is an optional step in the trajectory building
#' process before \code{preprocess_cds}.  After calculating dispersion for
#' a cell_data_set using the \code{calculate_gene_dispersion} function, the 
#' \code{select_genes} function allows the user to identify a set of genes
#' that will be used in downstream dimensionality reduction methods.  These
#' genes and their disperion and mean expression can be plotted using the 
#' \code{plot_gene_dispersion} function.
#'
#'
#' @param cds the cell_data_set upon which to perform this operation.
#' @param q the polynomial degree; default = 3.
#' @param id_tag the name of the feature data column corresponding to 
#' the unique id - typically ENSEMBL id; default = "id".
#' @param symbol_tag the name of the feature data column corresponding to 
#' the gene symbol; default = "gene_short_name".
#' @return an updated cell_data_set object with dispersion and mean expression saved
#' @export
#' 
calculate_gene_dispersion<-function(cds, q=3, id_tag="id", symbol_tag="gene_short_name"){
  m<-Matrix::rowMeans(counts(cds))
  # sd<-sqt(rowVars(as.matrix(counts(cds))))
  sd<-m3addon:::rowStdDev(exprs(cds))[1,]
  fdat<-fData(cds)
  cv<-sd/m*100
  df<-data.frame(log_dispersion=log(cv), log_mean=log(m))
  model <- lm(data = df, log_dispersion ~ log_mean + poly(log_mean, degree=q))
  prd <- data.frame(log_mean = df$log_mean)
  err<-suppressWarnings(predict(model, newdata= prd, se.fit = T))
  prd$lci <- err$fit - 1.96 * err$se.fit
  prd$fit <- err$fit
  prd$uci <- err$fit + 1.96 * err$se.fit
  prd$log_dispersion<-df$log_dispersion
  prd$id<-fdat[[id_tag]]
  prd$gene_short_name<-fdat[[symbol_tag]]
  cds@int_metadata$dispersion<-prd
  cds
}
