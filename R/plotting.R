#' Plot geneset scores
#'
#' @description Geneset scores are a score calculated for each single cell derived from \
#' more than one gene.  This function plots geneset scores using monocle3's 'plot_genes' function.
#' 
#' When using method 'totals', the sum of the size-factor corrected, log-normalized gene \
#' expression for a give set of genes is performed.  When using method 'corrected', single \
#' cell scores for a give gene set that have been "corrected" using 100X genes with similar \
#' expression levels.
 
#' @param cds Input cell_data_set object.
#' @param marker_set Vector of genes in the gene_metadata DataFrame (fData) that can be found in the column 'fData_col'
#' @param name Name given to the geneset
#' @param fData_col Character string denoting the gene_metadata DataFrame (fData) column that contains the genes in marker_set1.  Default = 'gene_short_name'
#' @return Plot
#' @importFrom Matrix colSums
#' @importFrom Matrix t
#' @export
#' @references Puram, S. V. et al. Single-Cell Transcriptomic Analysis of Primary and 
#' Metastatic Tumor Ecosystems in Head and Neck Cancer. Cell 171, 1611.e1–1611.e24 (2017).

plot_geneset<-function(cds, marker_set, name, fData_col="gene_short_name", method=c("totals", "corrected")){
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be one of 'totals' or 'corrected'")
  method <- match.arg(method)
  if(method=="totals") pData(cds)[[name]]<-estimate_score(cds, M2[[name]])
  if(method=="corrected") pData(cds)[[name]]<-estimate_corrected_score(cds, M2[[name]])
  nc<-nchar(name)
  if(nc>50){fontsize<-10}else{fontsize=14}
  switch(method, totals={loca="UW"}, 
         corrected={loca="Broad"})
  plot_cells(cds, color_cells_by = name, label_cell_groups = F, cell_size = 0.5)+ 
    theme(legend.position="top", legend.title = element_blank())+
    ggtitle(paste0(name, ": ", loca))+ 
    theme(plot.title = element_text(size = fontsize, face = "bold"), legend.text = element_text(size=9, angle = 90, vjust=0.5, hjust=0.3))+
    scale_color_gradientn(colors=c( "darkblue","skyblue", "white", "red", "darkred"))
}