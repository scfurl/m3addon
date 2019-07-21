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
    fData(cds)$exprs_bin = cut2(log(Matrix::rowMeans(normalized_counts(cds))), m=floor(nrow(fData(cds))/25))
  }
  cds
}


