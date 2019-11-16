remove_duplicated_rows<-function(cds, fdata_col="gene_short_name", unique_labels="id"){
  exprMat<-normalized_counts(cds)
  rn<-mcols(cds)[[fdata_col]]
  
  if (length(make.unique(rn)) != length(unique(rn))) {
    unique_rn<-rownames(exprMat)[!rn %in% rn[duplicated(rn)]]
    duplicated_rn<-rownames(exprMat)[rn %in% rn[duplicated(rn)]]
    duplicated_lab<-rn[rn %in% rn[duplicated(rn)]]
    warning("Non unique gene_labels found in data; these have been collapsed keeping the label that has the greatest variance of normalized data")
    duplicated<-exprMat[rn %in% rn[duplicated(rn)],]
    ts = data.table(rowV = rowVars(as.matrix(duplicated)), rn_D=duplicated_lab, rn_U=duplicated_rn)
    tokeep<-c(ts[ , .SD[which.max(rowV)], by = rn_D]$rn_U, unique_rn)
    rn<-rn[rownames(exprMat) %in% tokeep]
    cds[rownames(exprMat) %in% tokeep,]
  }
}


#' Pseudobulh=k
#'
#' @description This function creates pseudocells comprosed of a composite of data from n cells using 
#' an additive method.
#' @param cds Input cell_data_set object.
#' @param by column colData from which to group pseudobulking efforts
#' @param n number of cells to bulk
#' @return cds object
#' @importFrom monocle3 new_cell_data_set
#' @importFrom S4Vectors DataFrame
#' @importFrom Matrix Matrix
#' @export



pseudobulk<-function(cds, by, n=10){
  matlist<-lapply(levels(pData(cds)[[by]]), function(Cluster){
    datlist<-list()
    cds_ss<-cds[,pData(cds)[[by]] %in% Cluster]
    rand<-sample(1:length(colnames(cds_ss)), length(colnames(cds_ss)), replace=F)
    n <- ceiling(length(colnames(cds_ss))/n)
    pseudocelllist<-split(rand, 1:n)
    datlist[[1]]<-pseudocelllist
    datlist[[2]]<-lapply(pseudocelllist, function(pseudocell){
      Matrix::rowSums(normalized_counts(cds_ss[,pseudocell]))})
    names(datlist)<-c("indices", "pseudocells")
    return(datlist)
  })
  names(matlist)<-levels(pData(cds)[[by]])
  do.call(cbind, matlist[[1]]$pseudocells)
  #make new cds
  clusterlist<-lapply(matlist, function(alist) Matrix(do.call(cbind, alist$pseudocells), sparse=T))
  
  smat<-do.call(cbind, clusterlist)
  colnames(smat)<-unlist(sapply(1:length(matlist), function(num){
    pseudocellnum<-names(matlist[[num]]$indices)
    first_name<-paste(names(matlist)[num], pseudocellnum, sep="-")
    second_name<-sapply(matlist[[num]]$indices, paste, collapse=",")
    return(paste(first_name, second_name, sep="_"))
  }))
  
  combined_metadata_df<-data.frame(Pseudocell=sapply(strsplit(sapply(strsplit(colnames(smat), "-"), "[[",2), "_"), "[[", 1), row.names = colnames(smat))
  combined_metadata_df[[by]]=sapply(strsplit(colnames(smat), "-"), "[[",1)
  pd = DataFrame(combined_metadata_df)
  
  return(new_cell_data_set(smat*10, 
                 cell_metadata = pd, gene_metadata = fData(cds)))
  
}

#' ssGSEA
#'
#' @description This function runs ssGSEA as implemented here: https://gist.github.com/gaoce/39e0907146c752c127728ad74e123b33

#' @param X matrix. Rows are genes. Columns are samples. Row names are symbols.
#' @param gene_sets list. Each element is a string vector with gene symbols.
#' @param alpha numeric. Parameter for ssGSEA, the default is 0.25
#' @param scale logical. If True, normalize the scores by number of genes in the gene sets.
#' @param norm logical. If True, normalize the scores by the absolute difference between max and min values.
#' @param single logical. If True, use ssGSEA algorithm, otherwise use GSEA.
#'
#' @return matrix containing enrichment scroes. Rows are gene sets, columns are samples.
#' @importFrom matrixStats colRanks
#' @examples
#' # Create a fake matrix
#' m = 100
#' n = 100
#' set.seed(1)
#' X = matrix(rnorm(m*n), m, n)
#' # Assign 'gene symbols' to row names
#' rownames(X) = 1:m
#' # Create 3 gene sets
#' gene_sets = list(a = sample(m, 5), b = sample(m, 5), c = sample(m, 5))
#' system.time(assign('a', GSVA::gsva(X, gene_sets, method = 'ssgsea')))
#' system.time(assign('b', ssgsea(X, gene_sets, scale = F, norm = T)))
#' identical(a, b)
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  # Ranks for genes
  R = colRanks(X, preserveShape = T, ties.method = 'average')
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}



