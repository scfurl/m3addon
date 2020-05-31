#file<-"/Users/sfurla/OneDrive - Fred Hutchinson Cancer Research Center/computation/Analysis/ddata/thymus/data/SNF/subsampled_corrected.h5"
#' Imports a scanpy h5 object
#'
#' @description read in a Scanpy h5 object.
#'
#' @param file h5 scanpy object
#' @param raw data contained in the raw slot - See scanpy documentation for more info
#' @param rowname_col feature column for rownames (should contain unique values - suggest ENSID)
#' @param gene_short_name_col feature column for gene_short_name (suggest Symbol)
#' @param colors name or names of color palatte to import
#' @param exprs_bin Boolean whether to bin genes by mean expression
#' @param exprs_cuts Numeic indicating number of bins if using exprs_bin
#' @return Updated cell_data_set object
#' @import monocle3
#' @import Matrix
#' @import hdf5r
#' @importFrom S4Vectors DataFrame
#' @export
#' 
#' 

import_scanpy_h5<-function(file, raw=TRUE, rowname_col="GeneID", gene_short_name_col="GeneName", colors=NULL, PCA=T, UMAP=TRUE, var_features=T){
  if(!file.exists(file)){stop(paste0("Input file not found: ", file))}
  s <- H5File$new(file, mode = "r")
  if(raw){
    data = s[["raw"]][["X"]][["data"]][]
    indices = s[["raw/X/indices"]][] + 1
    indptr = s[["raw/X/indptr"]][]
    ncells<-length(indptr)-1
    obsdf<-DataFrame(rows=1:ncells)
  }else{
    data = s[["X/data"]][]
    indices = s[["X/indices"]][] + 1
    indptr = s[["X/indptr"]][]
    ncells<-length(indptr)-1
    obsdf<-DataFrame(rows=1:ncells)
  }
  for(name in s[["obs"]]$names){
    if(name=="__categories"){next}
    obsdf[[name]]<-s[["obs"]][[name]][]
  }
  for(name in names(s[["obs"]][["__categories"]])){
    obsdf[[name]]<-s[["obs"]][["__categories"]][[name]][][obsdf[[name]]+1]
  }
  
  nfeatures<-length(s[["raw/var/index"]][])
  gdf<-DataFrame(rows=1:nfeatures)
  for(name in s[["raw/var"]]$names){
    if(name=="__categories"){next}
    gdf[[name]]<-s[["raw/var"]][[name]][]
  }
  for(name in names(s[["raw/var"]][["__categories"]])){
    gdf[[name]]<-s[["raw/var"]][["__categories"]][[name]][][gdf[[name]]+1]
  }
  
  shape = c(nrow(obsdf), nrow(gdf))
  
  #umap<-read.csv(file.path(DATA_DIR, "thymus_annotated_matrix_files/umap.csv"), header=F)
  gbm = sparseMatrix(x = data, i = indices, p = indptr, dims = rev(shape))
  if("index" %in% colnames(obsdf)){
    cn<-obsdf$index
  }else{
    cn<-1:ncells
  }
  colnames(gbm) = cn
  rownames(obsdf)<-cn
  if(rowname_col %in% colnames(gdf)){
    rn<-gdf[[rowname_col]]
    if(any(duplicated(rn))){stop("Found duplicated values in suppliced 'rowname_col'")}
  }else{
    stop(paste0("Supplied 'rowname_col' not found in the gene metadata.  Found the following columns that may work for you:\n", 
                paste0(colnames(gdf), collapse="\n")))
  }
  rownames(gdf)<-gdf[[rowname_col]]
  rownames(gbm) = gdf[[rowname_col]]
  if(gene_short_name_col %in% colnames(gdf)){
    gdf$gene_short_name<-gdf$Symbol
  }else{
    warning(paste0("Supplied 'gene_short_name_col' not found in the gene metadata.  Found the following columns that may work for you:\n", 
                paste0(colnames(gdf), collapse="\n")))
  }

  suppressWarnings({
    cds = new_cell_data_set(gbm, cell_metadata = obsdf, 
                            gene_metadata = gdf)
  })
  
  if(!is.null(colors)){
    for( color in colors){
      cn<-paste0("Pallate_",color)
      cds@metadata[[cn]]<-s[["uns"]][[color]][]
    }
  }
  if(var_features) cds@metadata$var_features <- s[["var"]][[rowname_col]][]
  if(PCA) reducedDims(cds)<-SimpleList(PCA=t(s[["obsm"]][["X_pca"]]$read()))
  if(UMAP) reducedDims(cds)[["UMAP"]]=t(s[["obsm"]][["X_umap"]]$read())
  s$close()
  cds
  #plot_cells(cds, color_cells_by = "cell_type", label_cell_groups = F)+scale_color_manual(values=rev(cds@metadata$ct_col))
}
