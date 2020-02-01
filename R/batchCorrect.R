#' Performs mNN batch correction on PCA data from a cds
#'
#' @description 
#'
#' @param cds Input cell_data_set object.
#' @param batch_term variable in ColData (pData) of cds that informs the batch to be corrected
#' @param fast Whether to perform fastMNN or mnnCorrect (see scran package)
#' @param min_expr Numeric indicating expression threshold
#' @param exprs_bin Boolean whether to bin genes by mean expression
#' @param exprs_cuts Numeic indicating number of bins if using exprs_bin
#' @return Updated cell_data_set object
#' @importFrom batchelor fastMNN
#' @export
#' 
#' 

mnnCorrect_cds<-function(cds,  batch_term, fast=TRUE,
                        k=20, sigma=0.1, cos.norm.in=TRUE, cos.norm.out=TRUE,
                         svd.dim=0L, var.adj=TRUE, compute.angle=FALSE, subset.row=NULL, 
                         order=NULL, pc.approx=FALSE, irlba.args=list()){
  
  #get mnn args
  # args<-c(as.list(environment()))
  # args_to_pass<-c("k", "sigma", "cos.norm.in", "cos.norm.out",
  # "svd.dim", "var.adj", "compute.angle", "subset.row", 
  # "order", "pc.approx", "irlba.args")
  # mnn_args<-args[args_to_pass]
  # rm(args)
  #convert batch_term to factor if not already
  if(!class(pData(cds)[[batch_term]])=="factor"){
    pData(cds)[[batch_term]]<-factor(pData(cds)[[batch_term]])
  }
  
  # #get cell order
  # final_col_order<-colnames(cds)
  # 
  # #split PCA data into factor levels
  # mat_list<-list()
  # separated_colnames<-list()
  # i<-1
  # for(level in levels(pData(cds)[[batch_term]])){
  #   mat_list[[i]]<-t(reducedDims(cds)[["PCA"]][pData(cds)[[batch_term]] %in% level,])
  #   separated_colnames[[i]]<-colnames(mat_list[[i]])
  #   i<-i+1
  # }
  # i<-NULL
  # 
  
  # #combine matrices and mnn args
  # mnn_call<-c(mat_list, mnn_args)
  # 
  #call mnn
  
  if (fast){
    mnn_out<-fastMNN(t(reducedDims(cds)[["PCA"]]), batch=pData(cds)[[batch_term]])
  }else{
    mnn_out<-do.call(mnnCorrect, mnn_call)
  }
  
  mnn_rownames<-do.call(c, separated_colnames)
  corrected_mat<-mnn_out_corrected
  reducedDims(cds)[["PCA"]]<-corrected_mat[,match(mnn_rownames, original_colnames)]
}