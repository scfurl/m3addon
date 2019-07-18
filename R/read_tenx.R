#' @importFrom h5 h5file
read.cds.cellranger.h5.file = function(h5.file) {
  s<-h5file(h5.file)
  barcodes = readDataSet(s["matrix"]["barcodes"])
  gene_ids = readDataSet(s["matrix"]["features"]["id"])
  gene_names =readDataSet(s["matrix"]["features"]["name"])
  data = as.double(readDataSet(s["matrix"]["data"]))
  indices = as.integer(readDataSet(s["matrix"]["indices"]))
  indptr = as.integer(readDataSet(s["matrix"]["indptr"]))


  gbm = new(
    "dgCMatrix",
    x = data, i = indices, p = indptr,
    Dim = c(length(gene_ids), length(barcodes))
  )

  pData.df = data.frame(
    barcode = barcodes,
    stringsAsFactors = F
  )

  fData.df = data.frame(
    id = gene_ids,
    gene_short_name = gene_names,
    stringsAsFactors = F
  )

  rownames(pData.df) = pData.df$barcode
  rownames(fData.df) = fData.df$id

  rownames(gbm) = rownames(fData.df)
  colnames(gbm) = rownames(pData.df)

  suppressWarnings({res = new_cell_data_set(
    gbm,
    cell_metadata = pData.df,
    gene_metadata = fData.df
  )})

  return(res)
}

#' Make monocle3 cell_data_set using cellranger h5
#' @description This function reads a vector of cellranger "out" folders for .h5 files.  This function
#' will create a cds of the data using cellranger thresholds for filtering out droplets that do not contain cells.
#' Optionally, this function can return the unfiltered data as well for further manipulation (ie adjustment of 
#' droplet thresholds) if desired.
#'
#' @param folders A vector of cellranger folders (full folder name is ideal).  These must contain an "outs" subfolders with
#' two files; 1) raw_feature_bc_matrix.h5, and 2) filtered_feature_bc_matrix.h5.  For aggregated samples
#' run with cellrangers 'aggregate' function, this folder should also contain the aggregation.csv file.
#' @param samplenames An optional vector that corresponds to the names you would like to give to
#'each element in filelist.  This defaults to the the basename of the filelist element i.e.
# the folder /home/user/cellranger/runs/Samp1-SA-GA_A1 would be given the name: "Samp1-SA-GA_A1"
#' @param unfiltered This parameter, if true, returns a list containing 1) the filtered cds,
#'2) a list of unfiltered cds for each sample in 'folders'.
#' @param empty.droplet.threshold minimum number of umis per droplet to include in unfiltered output (default 15)
#' @param expressed_genes If true, this option removes genes that do not have expression in at least a 
#' minumum number of cells (cell_min parameter)
#' @param cell_min Minimum number of cells that a gene needs to be expressed for the gene to be included in the 
#' @param aggregated will read a cellranger aggregate folder, this is currently only supported for 1 folder.  The 'aggregation.csv'
#' file must be included in the 'outs' folder.
#' @return a cell_data_set object or a list of items if unfiltered data is returned (see unfiltered)
#' @importFrom h5 h5file
#' @export
load_cellranger_data_h5<-function(folders, 
                                  samplenames=NULL, 
                                  unfiltered=F, 
                                  empty.droplet.threshold=15, 
                                  expressed_genes=TRUE,
                                  cell_min=5,
                                  aggregated=F
                                  ){
  #This function reads a vector of cellranger folder "out" folders for .h5 files.  Optionally can
  #return the unfiltered data as well.  
  if(aggregated){
    if(length(folders)>1) stop("Only one aggregated folder supported")
    message(paste0("Reading aggregated data for: ", folders[1]))
    if(!file.exists(file.path(folders, "outs", "aggregation.csv"))) stop("aggregation.csv file must be present in 'outs'")
    cds = read.cds.cellranger.h5.file(
      file.path(folders[1], "outs", "raw_feature_bc_matrix.h5"))
    agg<-read.csv(file.path(folders[1], "outs", "aggregation.csv"))
    message(paste0("Finding threshold for aggregegated data: ", folders[1]))
    pData(cds)$n.umi = Matrix::colSums(exprs(cds))
    h5.file<-file.path(folders[1], "outs", "filtered_feature_bc_matrix.h5")
    s<-h5file(h5.file)
    umi.count.thresholds.from.cellranger = pData(cds)$n.umi[order(-pData(cds)$n.umi)][readDataSet(s["matrix"]["shape"])[2]]
    umi.count.thresholds = umi.count.thresholds.from.cellranger
    cds = cds[,pData(cds)$n.umi >= umi.count.thresholds                                                                ]
    pData(cds)$sample_no<-sapply(strsplit(rownames(pData(cds)), "-"), "[[", 2)
    agg$sample_no<-1:nrow(agg)
    if(expressed_genes){
      cds<-detect_genes(cds)
      cds<-cds[fData(cds)$num_cells_expressed>cell_min,]
    }
    cds<-add_meta_data_cds(cds=cds, meta=agg, cds_col = "sample_no", meta_col = "sample_no")
    return(cds)
  }
  if(is.null(samplenames)){
    sample.ids<-folders
    names(sample.ids)<-sapply(folders, basename)
  }else{
    sample.ids<-folders
    names(sample.ids)<-samplenames
  }

  umi.count.thresholds.from.cellranger = list()
  unfiltered.cds.list<-list()
  umi.count.thresholds<-list()


  #read data
  for(sample.id in sample.ids){
    message(paste0("Reading data for: ", sample.id))
    unfiltered.cds.list[[sample.id]] = read.cds.cellranger.h5.file(
      file.path(sample.id, "outs", "raw_feature_bc_matrix.h5"))
  }

  #Get threshes from cellranger
  for (sample.id in sample.ids) {
    message(paste0("Finding threshold for: ", sample.id))
    pData(unfiltered.cds.list[[sample.id]])$n.umi =
      Matrix::colSums(exprs(unfiltered.cds.list[[sample.id]]))

    # umi.cdf = ecdf(pData(unfiltered.cds.list[[sample.id]])$n.umi)
    #
    # pData(unfiltered.cds.list[[sample.id]])$umi.quantile =
    #   umi.cdf(pData(unfiltered.cds.list[[sample.id]])$n.umi)

    h5.file<-file.path(sample.id, "outs", "filtered_feature_bc_matrix.h5")
    s<-h5file(h5.file)
    umi.count.thresholds.from.cellranger[[sample.id]] = pData(unfiltered.cds.list[[sample.id]])$n.umi[order(-pData(unfiltered.cds.list[[sample.id]])$n.umi)][readDataSet(s["matrix"]["shape"])[2]]

    umi.count.thresholds[[sample.id]] = umi.count.thresholds.from.cellranger[[sample.id]]
  }

  
  filtered.cds.list = list()
  for (sample.id in sample.ids) {
    filtered.cds.list[[sample.id]] = unfiltered.cds.list[[sample.id]][,
                                                                      pData(unfiltered.cds.list[[sample.id]])$n.umi >=
                                                                        umi.count.thresholds[[sample.id]]
                                                                      ]
  }

  #sapply(sample.ids, function(x) ncol(filtered.cds.list[[x]]))
  #rm(unfiltered.cds.list)
  fdat_rownames<-lapply(filtered.cds.list, function(cds) rownames(fData(cds)))
  if(!all.identical(fdat_rownames))stop("Not all genes are the same across samples")
  common.fData = fData(filtered.cds.list[[sample.ids[1]]])
  names(unfiltered.cds.list)<-names(sample.ids)
  names(filtered.cds.list)<-names(sample.ids)
  new.pData = list()
  for (sample.id in names(sample.ids)) {
    new.pData[[sample.id]] = pData(filtered.cds.list[[sample.id]])

    new.pData[[sample.id]]$sample = sample.id

    new.pData[[sample.id]]$cell = paste(
      sample.id, new.pData[[sample.id]]$barcode, sep = ".")

    rownames(new.pData[[sample.id]]) = new.pData[[sample.id]]$cell
  }

  new.exprs = list()
  for (sample.id in names(sample.ids)) {
    new.exprs[[sample.id]] = exprs(filtered.cds.list[[sample.id]])
    colnames(new.exprs[[sample.id]]) = new.pData[[sample.id]]$cell
  }

  combined.pData = do.call(rbind, new.pData)
  combined.pData = combined.pData[, c("barcode", "n.umi",  "sample")]
  rownames(combined.pData) = paste0(combined.pData$sample, ".", combined.pData$barcode)
  #head(combined.pData, 2)


  combined.exprs = do.call(cbind, new.exprs)


  #colnames(combined.exprs)==rownames(combined.pData)

  cds = new_cell_data_set(
    combined.exprs,
    cell_metadata =  combined.pData,
    gene_metadata = common.fData
  )
  if(unfiltered) return(list(cds=cds, threshes=umi.count.thresholds, unfiltered=unfiltered.cds.list))
  if(expressed_genes){
    cds<-detect_genes(cds)
    cds<-cds[fData(cds)$num_cells_expressed>cell_min,]
  }
  cds
}

#' Add colData (aka pData) to a cell_data_set
#' @description This function will take a cds and add sample levels phenodata to all the cells from that 
#sample in the cds.  The function will match sample names from the cds_anchor parameter to 
#the sample names from the df_anchor parameter and add the remaining phenodata from the df 
#to the cds.
#'
#' @param folders A vector of cellranger folders (full folder name is ideal).  These must contain an "outs" subfolders with
#' two files; 1) raw_feature_bc_matrix.h5, and 2) filtered_feature_bc_matrix.h5.  For aggregated samples
#' run with cellrangers 'aggregate' function, this folder should also contain the aggregation.csv file.
#' @param samplenames An optional vector that corresponds to the names you would like to give to
#'each element in filelist.  This defaults to the the basename of the filelist element i.e.
# the folder /home/user/cellranger/runs/Samp1-SA-GA_A1 would be given the name: "Samp1-SA-GA_A1"
#' @param unfiltered This parameter, if true, returns a list containing 1) the filtered cds,
#'2) a list of unfiltered cds for each sample in 'folders'.
#' @param empty.droplet.threshold minimum number of umis per droplet to include in unfiltered output (default 15)
#' @param expressed_genes If true, this option removes genes that do not have expression in at least a 
#' minumum number of cells (cell_min parameter)
#' @param cell_min Minimum number of cells that a gene needs to be expressed for the gene to be included in the 
#' @param aggregated will read a cellranger aggregate folder, this is currently only supported for 1 folder.  The 'aggregation.csv'
#' file must be included in the 'outs' folder.
#' @return a cell_data_set object or a list of items if unfiltered data is returned (see unfiltered)
#' @export
add_meta_data_cds<-function(cds, meta, cds_col="sample", meta_col="sample"){

  cds_col<-as.character(cds_col)
  meta_col<-as.character(meta_col)
  if(class(cds)!="cell_data_set"){stop("input 'cds' is not a monocle3 cell data set")}
  if(length(which(colnames(pData(cds)) %in% cds_col))!=1){stop("cds_col doesnt match only once")}
  if(length(which(colnames(meta) %in% meta_col))!=1){stop("meta_col doesnt match only once")}
  cds_ind<-match(cds_col, colnames(pData(cds)))
  to_add<-meta[match(as.character(pData(cds)[,cds_ind]), as.character(meta[,meta_col])),]
  meta_ind<-match(meta_col, colnames(meta))
  add<-to_add[,-meta_ind]
  rownames(add)<-rownames(pData(cds))
  pData(cds)<-cbind(pData(cds), add)
  return(cds)
}


all.identical <- function(l) all(mapply(identical, head(l, 1), tail(l, -1)))

