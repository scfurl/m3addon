
#' @import hdf5r
#' @importFrom Matrix sparseMatrix
read.cds.cellranger.h5.file = function(h5.file) {
  if(!file.exists(h5.file)){stop(paste0("File ", h5.file," not found"))}
  #s<-h5file(h5.file, mode="r")
  s <- H5File$new(h5.file, mode="r+")

  #s$ls(recursive=T)

  barcodes = s[["matrix/barcodes"]][]
  gene_ids = s[["matrix/features/id"]][]
  gene_names =s[["matrix/features/name"]][]
  featuretype = s[["matrix/features/feature_type"]][]
  data = s[["matrix/data"]][]
  indices = s[["matrix/indices"]][]+1
  indptr = s[["matrix/indptr"]][]
  shape = s[["matrix/shape"]][]
  #browser()
  h5close(s)
  # gbm = new(
  #   "dgCMatrix",
  #   x = data, i = indices, p = indptr,
  #   Dim = shape)
  if(length(levels(factor(featuretype)))>1){
    warning("Warning - Multiple feature types found")
  }
  
  gbm = sparseMatrix(x = data, i = indices, p = indptr,
    dims = shape)

  
  pData.df = data.frame(
    barcode = barcodes,
    stringsAsFactors = F
  )

  fData.df = data.frame(
    id = gene_ids,
    gene_short_name = gene_names,
    feature_type = featuretype,
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



#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix readMM
#' @importFrom data.table fread
read.cds.starsolo.file = function(folder) {
  if(!file.exists(folder)){stop(paste0("File ", folder," not found"))}
  files<-c("features.tsv", "barcodes.tsv", "matrix.mtx")
  
  barcodes = fread(file.path(folder, files[2]), header = F)$V1
  feats<- fread(file.path(folder, files[1]), header=F)
  gene_ids = feats$V1
  gene_names = feats$V2
  gbm = readMM(file.path(folder, files[3]))
  
  
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
#' @param folders A vector of cellranger folders (full folder name is ideal).  These must contain an "outs" subfolders with
#' two files; 1) raw_feature_bc_matrix.h5, and 2) filtered_feature_bc_matrix.h5.  For aggregated samples
#' run with cellrangers 'aggregate' function, this folder should also contain the aggregation.csv file.
#' @param samplenames An optional vector that corresponds to the names you would like to give to
#'each element in filelist.  This defaults to the the basename of the filelist element i.e.
# the folder /home/user/cellranger/runs/Samp1-SA-GA_A1 would be given the name: "Samp1-SA-GA_A1"
#' @param files A vector of filtered h5 cellranger files.  Not compatible with aggregated mode.
#' @param unfilt_files A vector of unfiltered h5 cellranger files.  Not compatible with aggregated mode.  Will be ignored if unfiltered is FALSE.
#' @param unfiltered This parameter, if true, returns a list containing 1) the filtered cds,
#'2) a list of unfiltered cds for each sample in 'folders'.
#' @param empty.droplet.threshold minimum number of umis per droplet to include in unfiltered output (default 15)
#' @param expressed_genes If true, this option removes genes that do not have expression in at least a 
#' minumum number of cells (cell_min parameter)
#' @param cell_min Minimum number of cells that a gene needs to be expressed for the gene to be included in the 
#' @param aggregated will read a cellranger aggregate folder, this is currently only supported for 1 folder.  The 'aggregation.csv'
#' file must be included in the 'outs' folder.
#' @return a cell_data_set object or a list of items if unfiltered data is returned (see unfiltered)
#' @importFrom Matrix colSums
#' @export
load_cellranger_data_h5<-function(folders=NULL, 
                                  samplenames=NULL, 
                                  unfiltered=F,
                                  files=NULL, unfilt_files=NULL,
                                  empty.droplet.threshold=15, 
                                  expressed_genes=TRUE,
                                  cell_min=1,
                                  aggregated=F, chemistry = "SC3Pv3", atac_feature="peaks"
                                  ){
  #This function reads a vector of cellranger folder "out" folders for .h5 files.  Optionally can
  #return the unfiltered data as well.  
  if(!is.null(files) & !is.null(folders)){stop("Run either specifying folders or files")}
  if(!is.null(files) & aggregated){stop("File mode not compatible with aggregated")}
  if(!is.null(files) & is.null(unfilt_files) & unfiltered){stop("Running unfiltered files, but filenames (unfilt_files) not provided")}
  if(!is.null(files)){filemode<-T}else{filemode<-F}
  if(!chemistry %in% c("threeprime", "fiveprime", "SC3Pv1", "SC3Pv2", "SC3Pv3", "SC5P-PE", "SC5P-R2", "ATAC")){stop("Chemistry not found")}
  if(chemistry %in% c("threeprime", "fiveprime", "SC3Pv1", "SC3Pv2", "SC3Pv3", "SC5P-PE", "SC5P-R2")){
    if(aggregated){
      #aggregate filtered option:
      if(length(folders)>1) stop("Only one aggregated folder supported")
      filt_file<-file.path(folders[1], "outs", "filtered_feature_bc_matrix.h5")
      unfilt_file<-file.path(folders[1], "outs", "raw_feature_bc_matrix.h5")
      agg_file<-file.path(folders[1], "outs", "aggregation.csv")
      message(paste0("Reading aggregated (filtered) data for: ", filt_file))
      if(!file.exists(agg_file)) stop("aggregation.csv file must be present in 'outs'")
      cds = read.cds.cellranger.h5.file(filt_file)
      agg<-read.csv(agg_file)
      pData(cds)$n.umi = Matrix::colSums(exprs(cds))
      pData(cds)$sample_no<-sapply(strsplit(rownames(pData(cds)), "-"), "[[", 2)
      agg$sample_no<-1:nrow(agg)
      if(expressed_genes){
        cds<-detect_genes(cds, exprs_bin = F)
        cds<-cds[fData(cds)$num_cells_expressed>cell_min,]
      }
      cds<-add_meta_data_cds(cds=cds, meta=agg, cds_col = "sample_no", meta_col = "sample_no")
      if(aggregated & !unfiltered){return(cds)}
      
      #aggregate unfiltered option:
      message(paste0("Reading aggregated (unfiltered) data for: ", filt_file))
      cds_unfilt = read.cds.cellranger.h5.file(unfilt_file)
      agg<-read.csv(agg_file)
      pData(cds_unfilt)$n.umi = Matrix::colSums(exprs(cds_unfilt))
      pData(cds_unfilt)$sample_no<-sapply(strsplit(rownames(pData(cds_unfilt)), "-"), "[[", 2)
      agg$sample_no<-1:nrow(agg)
      if(expressed_genes){
        cds_unfilt<-detect_genes(cds_unfilt, exprs_bin = F)
        cds_unfilt<-cds_unfilt[fData(cds_unfilt)$num_cells_expressed>cell_min,]
      }
      cds_unfilt<-add_meta_data_cds(cds=cds_unfilt, meta=agg, cds_col = "sample_no", meta_col = "sample_no")
      if(aggregated & unfiltered){return(list(unfiltered_cds=cds_unfilt, filtered_cds=cds))}
    }
  
    #browser()
    #multiple files (no_agg); unfiltered option
    if(filemode){
      filt_folders<-files
    }else{
      filt_folders<-file.path(folders, "outs", "filtered_feature_bc_matrix.h5")
    }
    if(is.null(samplenames)){
      sample.ids<-filt_folders
      names(sample.ids)<-sapply(filt_folders, basename)
    }else{
      sample.ids<-filt_folders
      names(sample.ids)<-samplenames
    }
    
    ######filt
    #read filtered data
    filtered.cds.list<-list()
    
    for(sample.id in sample.ids){
      message(paste0("Reading (filtered) data for: ", sample.id))
      filtered.cds.list[[sample.id]] = read.cds.cellranger.h5.file(sample.id)
      pData(filtered.cds.list[[sample.id]])$n.umi<-colSums(exprs(filtered.cds.list[[sample.id]]))
    }
    
    #add_and checkrownames
    if(length(filtered.cds.list)>1){
      fdat_rownames<-lapply(filtered.cds.list, function(cds) rownames(fData(cds)))
      if(!all.identical(fdat_rownames))stop("Not all genes are the same across samples")
    }
    #make fData
    common.fData = fData(filtered.cds.list[[sample.ids[1]]])
    names(filtered.cds.list)<-names(sample.ids)
    #make pData
    new.pData = list()
    for (sample.id in names(sample.ids)) {
      new.pData[[sample.id]] = pData(filtered.cds.list[[sample.id]])
      
      new.pData[[sample.id]]$sample = sample.id
      
      new.pData[[sample.id]]$cell = paste(
        sample.id, new.pData[[sample.id]]$barcode, sep = ".")
      
      rownames(new.pData[[sample.id]]) = new.pData[[sample.id]]$cell
    }
    #make exprsData
    new.exprs = list()
    for (sample.id in names(sample.ids)) {
      new.exprs[[sample.id]] = exprs(filtered.cds.list[[sample.id]])
      colnames(new.exprs[[sample.id]]) = new.pData[[sample.id]]$cell
    }
    
    #combine
    combined.pData = do.call(rbind, new.pData)
    combined.pData = combined.pData[, c("barcode", "n.umi",  "sample")]
    rownames(combined.pData) = paste0(combined.pData$sample, ".", combined.pData$barcode)
    combined.exprs = do.call(cbind, new.exprs)
    cds = new_cell_data_set(
      combined.exprs,
      cell_metadata =  combined.pData,
      gene_metadata = common.fData
    )
    if(expressed_genes){
      cds<-detect_genes(cds, exprs_bin = F)
      cds<-cds[fData(cds)$num_cells_expressed>cell_min,]
    }
    if(!unfiltered){return(cds)}
  
    #####unfilt
    #read unfiltered data
    unfiltered.cds.list<-list()
    if(filemode){
      unfilt_folders<-unfilt_files
    }else{
      unfilt_folders<-file.path(folders, "outs", "raw_feature_bc_matrix.h5")
    }
    if(is.null(samplenames)){
      sample.ids<-unfilt_folders
      names(sample.ids)<-sapply(unfilt_folders, basename)
    }else{
      sample.ids<-unfilt_folders
      names(sample.ids)<-samplenames
    }
    for(sample.id in sample.ids){
        message(paste0("Reading (unfiltered) data for: ", sample.id))
        unfiltered.cds.list[[sample.id]] = read.cds.cellranger.h5.file(
          file.path(sample.id))
        pData(unfiltered.cds.list[[sample.id]])$n.umi<-colSums(exprs(unfiltered.cds.list[[sample.id]]))
        #filter out based on empty droplet threshold
        unfiltered.cds.list[[sample.id]]<-unfiltered.cds.list[[sample.id]][,unfiltered.cds.list[[sample.id]]$n.umi>empty.droplet.threshold]
    }
    
    #add_and checkrownames
    if(length(filtered.cds.list)>1){
      fdat_rownames<-lapply(filtered.cds.list, function(cds) rownames(fData(cds)))
      if(!all.identical(fdat_rownames))stop("Not all genes are the same across samples")
    }
    #make fData
    common.fData = fData(unfiltered.cds.list[[sample.ids[1]]])
    names(unfiltered.cds.list)<-names(sample.ids)
    #make pData
    new.pData = list()
    for (sample.id in names(sample.ids)) {
      new.pData[[sample.id]] = pData(unfiltered.cds.list[[sample.id]])
      
      new.pData[[sample.id]]$sample = sample.id
      
      new.pData[[sample.id]]$cell = paste(
        sample.id, new.pData[[sample.id]]$barcode, sep = ".")
      
      rownames(new.pData[[sample.id]]) = new.pData[[sample.id]]$cell
    }
    #make exprsData
    new.exprs = list()
    for (sample.id in names(sample.ids)) {
      new.exprs[[sample.id]] = exprs(unfiltered.cds.list[[sample.id]])
      colnames(new.exprs[[sample.id]]) = new.pData[[sample.id]]$cell
    }
    
    #combine
    combined.pData = do.call(rbind, new.pData)
    combined.pData = combined.pData[, c("barcode", "n.umi",  "sample")]
    rownames(combined.pData) = paste0(combined.pData$sample, ".", combined.pData$barcode)
    combined.exprs = do.call(cbind, new.exprs)
    cds_unfilt = new_cell_data_set(
      combined.exprs,
      cell_metadata =  combined.pData,
      gene_metadata = common.fData
    )
    if(expressed_genes){
      cds_unfilt<-detect_genes(cds_unfilt, exprs_bin = F)
      cds_unfilt<-cds_unfilt[fData(cds_unfilt)$num_cells_expressed>cell_min,]
    }
    if(unfiltered){return(list(unfiltered_cds=cds_unfilt, filtered_cds=cds))}
  }
  if(chemistry %in% "ATAC"){
    if(aggregated){stop("Not currently supported")}
    if(unfiltered){stop("Not currently supported")}
    #browser()
    #multiple files (no_agg); unfiltered option
    if(is.null(samplenames)){
      sample.ids<-folders
      names(sample.ids)<-sapply(folders, basename)
    }else{
      sample.ids<-folders
      names(sample.ids)<-samplenames
    }
    if(!atac_feature %in% c("peaks", "tfs"))stop("Currently only peak and tf data supported")
    if(atac_feature=="peaks"){h5file<-"filtered_peak_bc_matrix.h5"}
    if(atac_feature=="tfs"){h5file<-"filtered_tf_bc_matrix.h5"}
    ######filt
    #read filtered data
    filtered.cds.list<-list()
    for(sample.id in sample.ids){
      message(paste0("Reading (filtered) data for: ", sample.id))
      filtered.cds.list[[sample.id]] = read.cds.cellranger.h5.file(
        file.path(sample.id, "outs", h5file))
      pData(filtered.cds.list[[sample.id]])$n.umi<-colSums(exprs(filtered.cds.list[[sample.id]]))
    }
    
    #add_and checkrownames
    if(length(filtered.cds.list)>1){
      fdat_rownames<-lapply(filtered.cds.list, function(cds) rownames(fData(cds)))
      if(!all.identical(fdat_rownames))stop("Not all genes are the same across samples")
    }
    #make fData
    common.fData = fData(filtered.cds.list[[sample.ids[1]]])
    names(filtered.cds.list)<-names(sample.ids)
    #make pData
    new.pData = list()
    for (sample.id in names(sample.ids)) {
      new.pData[[sample.id]] = pData(filtered.cds.list[[sample.id]])
      
      new.pData[[sample.id]]$sample = sample.id
      
      new.pData[[sample.id]]$cell = paste(
        sample.id, new.pData[[sample.id]]$barcode, sep = ".")
      
      rownames(new.pData[[sample.id]]) = new.pData[[sample.id]]$cell
    }
    #make exprsData
    new.exprs = list()
    for (sample.id in names(sample.ids)) {
      new.exprs[[sample.id]] = exprs(filtered.cds.list[[sample.id]])
      colnames(new.exprs[[sample.id]]) = new.pData[[sample.id]]$cell
    }
    
    #combine
    combined.pData = do.call(rbind, new.pData)
    combined.pData = combined.pData[, c("barcode", "n.umi",  "sample")]
    rownames(combined.pData) = paste0(combined.pData$sample, ".", combined.pData$barcode)
    combined.exprs = do.call(cbind, new.exprs)
    cds = new_cell_data_set(
      combined.exprs,
      cell_metadata =  combined.pData,
      gene_metadata = common.fData
    )
    if(expressed_genes){
      cds<-detect_genes(cds, exprs_bin = F)
      cds<-cds[fData(cds)$num_cells_expressed>cell_min,]
    }
    if(!unfiltered){return(cds)}
    
    #####unfilt
    #read unfiltered data
    unfiltered.cds.list<-list()
    for(sample.id in sample.ids){
      message(paste0("Reading (unfiltered) data for: ", sample.id))
      unfiltered.cds.list[[sample.id]] = read.cds.cellranger.h5.file(
        file.path(sample.id, "outs", "raw_feature_bc_matrix.h5"))
      pData(unfiltered.cds.list[[sample.id]])$n.umi<-colSums(exprs(unfiltered.cds.list[[sample.id]]))
    }
    
    #add_and checkrownames
    if(length(filtered.cds.list)>1){
      fdat_rownames<-lapply(filtered.cds.list, function(cds) rownames(fData(cds)))
      if(!all.identical(fdat_rownames))stop("Not all genes are the same across samples")
    }
    #make fData
    common.fData = fData(unfiltered.cds.list[[sample.ids[1]]])
    names(unfiltered.cds.list)<-names(sample.ids)
    #make pData
    new.pData = list()
    for (sample.id in names(sample.ids)) {
      new.pData[[sample.id]] = pData(unfiltered.cds.list[[sample.id]])
      
      new.pData[[sample.id]]$sample = sample.id
      
      new.pData[[sample.id]]$cell = paste(
        sample.id, new.pData[[sample.id]]$barcode, sep = ".")
      
      rownames(new.pData[[sample.id]]) = new.pData[[sample.id]]$cell
    }
    #make exprsData
    new.exprs = list()
    for (sample.id in names(sample.ids)) {
      new.exprs[[sample.id]] = exprs(unfiltered.cds.list[[sample.id]])
      colnames(new.exprs[[sample.id]]) = new.pData[[sample.id]]$cell
    }
    
    #combine
    combined.pData = do.call(rbind, new.pData)
    combined.pData = combined.pData[, c("barcode", "n.umi",  "sample")]
    rownames(combined.pData) = paste0(combined.pData$sample, ".", combined.pData$barcode)
    combined.exprs = do.call(cbind, new.exprs)
    cds_unfilt = new_cell_data_set(
      combined.exprs,
      cell_metadata =  combined.pData,
      gene_metadata = common.fData
    )
    if(expressed_genes){
      cds_unfilt<-detect_genes(cds_unfilt, exprs_bin = F)
      cds_unfilt<-cds_unfilt[fData(cds_unfilt)$num_cells_expressed>cell_min,]
    }
    if(unfiltered){return(list(unfiltered_cds=cds_unfilt, filtered_cds=cds))}
  }
}


#' Make monocle3 cell_data_set using STARsolo
#' @description This function reads a vector of STARsolo "out" folders for .mtx and .tsv files.  This function
#' will create a cds of the data using cellranger thresholds for filtering out droplets that do not contain cells.
#' Optionally, this function can return the unfiltered data as well for further manipulation (ie adjustment of 
#' droplet thresholds) if desired.
#' @param folders A vector of STARsolo folders (full folder name is ideal).  These must contain an "outs" subfolders with
#' three files: 1) features.tsv, 2) barcodes.tsv, and 3) matrix.mtx.
#' @param samplenames An optional vector that corresponds to the names you would like to give to
#'each element in filelist.
#' @param count_method (required) - A string denoting STARsolo output counting method.  (one of either: Gene, GeneFull, SJ, Velocyto)
#' @param unfiltered This parameter, if true, returns a list containing 1) the filtered cds,
#'2) a list of unfiltered cds for each sample in 'folders'.
#' @param empty.droplet.threshold minimum number of umis per droplet to include in unfiltered output (default 15)
#' @param expressed_genes If true, this option removes genes that do not have expression in at least a 
#' minumum number of cells (cell_min parameter)
#' @param cell_min Minimum number of cells that a gene needs to be expressed for the gene to be included in the 
#' @return a cell_data_set object or a list of items if unfiltered data is returned (see unfiltered)
#' @importFrom Matrix colSums
#' @export
load_STARsolo_data<-function(folders, 
                                  samplenames=NULL, 
                                  unfiltered=F, 
                                  empty.droplet.threshold=15, 
                                  expressed_genes=TRUE, count_method="Gene",
                                  cell_min=1){
    

    if(is.null(samplenames)){
      sample.ids<-folders
      names(sample.ids)<-sapply(folders, basename)
    }else{
      sample.ids<-folders
      names(sample.ids)<-samplenames
    }
    
    ######filt
    #read filtered data
    filtered.cds.list<-list()
    for(sample.id in sample.ids){
      message(paste0("Reading (filtered) data for: ", sample.id))
      filtered.cds.list[[sample.id]] = read.cds.starsolo.file(
        file.path(sample.id, count_method, "filtered"))
      pData(filtered.cds.list[[sample.id]])$n.umi<-colSums(exprs(filtered.cds.list[[sample.id]]))
    }
    
    #add_and checkrownames
    if(length(filtered.cds.list)>1){
      fdat_rownames<-lapply(filtered.cds.list, function(cds) rownames(fData(cds)))
      if(!all.identical(fdat_rownames))stop("Not all genes are the same across samples")
    }
    #make fData
    common.fData = fData(filtered.cds.list[[sample.ids[1]]])
    names(filtered.cds.list)<-names(sample.ids)
    #make pData
    new.pData = list()
    for (sample.id in names(sample.ids)) {
      new.pData[[sample.id]] = pData(filtered.cds.list[[sample.id]])
      
      new.pData[[sample.id]]$sample = sample.id
      
      new.pData[[sample.id]]$cell = paste(
        sample.id, new.pData[[sample.id]]$barcode, sep = ".")
      
      rownames(new.pData[[sample.id]]) = new.pData[[sample.id]]$cell
    }
    #make exprsData
    new.exprs = list()
    for (sample.id in names(sample.ids)) {
      new.exprs[[sample.id]] = exprs(filtered.cds.list[[sample.id]])
      colnames(new.exprs[[sample.id]]) = new.pData[[sample.id]]$cell
    }
    
    #combine
    combined.pData = do.call(rbind, new.pData)
    combined.pData = combined.pData[, c("barcode", "n.umi",  "sample")]
    rownames(combined.pData) = paste0(combined.pData$sample, ".", combined.pData$barcode)
    combined.exprs = do.call(cbind, new.exprs)
    cds = new_cell_data_set(
      combined.exprs,
      cell_metadata =  combined.pData,
      gene_metadata = common.fData
    )
    if(expressed_genes){
      cds<-detect_genes(cds, exprs_bin = F)
      cds<-cds[fData(cds)$num_cells_expressed>cell_min,]
    }
    if(!unfiltered){return(cds)}
    
    #####unfilt
    #read unfiltered data
    unfiltered.cds.list<-list()
    for(sample.id in sample.ids){
      message(paste0("Reading (unfiltered) data for: ", sample.id))
      unfiltered.cds.list[[sample.id]] = read.cds.starsolo.file(
        file.path(sample.id, count_method, "raw"))
      pData(unfiltered.cds.list[[sample.id]])$n.umi<-colSums(exprs(unfiltered.cds.list[[sample.id]]))
    }
    
    #add_and checkrownames
    if(length(filtered.cds.list)>1){
      fdat_rownames<-lapply(filtered.cds.list, function(cds) rownames(fData(cds)))
      if(!all.identical(fdat_rownames))stop("Not all genes are the same across samples")
    }
    #make fData
    common.fData = fData(unfiltered.cds.list[[sample.ids[1]]])
    names(unfiltered.cds.list)<-names(sample.ids)
    #make pData
    new.pData = list()
    for (sample.id in names(sample.ids)) {
      new.pData[[sample.id]] = pData(unfiltered.cds.list[[sample.id]])
      
      new.pData[[sample.id]]$sample = sample.id
      
      new.pData[[sample.id]]$cell = paste(
        sample.id, new.pData[[sample.id]]$barcode, sep = ".")
      
      rownames(new.pData[[sample.id]]) = new.pData[[sample.id]]$cell
    }
    #make exprsData
    new.exprs = list()
    for (sample.id in names(sample.ids)) {
      new.exprs[[sample.id]] = exprs(unfiltered.cds.list[[sample.id]])
      colnames(new.exprs[[sample.id]]) = new.pData[[sample.id]]$cell
    }
    
    #combine
    combined.pData = do.call(rbind, new.pData)
    combined.pData = combined.pData[, c("barcode", "n.umi",  "sample")]
    rownames(combined.pData) = paste0(combined.pData$sample, ".", combined.pData$barcode)
    combined.exprs = do.call(cbind, new.exprs)
    cds_unfilt = new_cell_data_set(
      combined.exprs,
      cell_metadata =  combined.pData,
      gene_metadata = common.fData
    )
    if(expressed_genes){
      cds_unfilt<-detect_genes(cds_unfilt, exprs_bin = F)
      cds_unfilt<-cds_unfilt[fData(cds_unfilt)$num_cells_expressed>cell_min,]
    }
    if(unfiltered){return(list(unfiltered_cds=cds_unfilt, filtered_cds=cds))}
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

