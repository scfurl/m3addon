#' @title Doublet Detection
#' @description See Gayoso, Adam, & Shor, Jonathan. (2018, July 17). 
#' DoubletDetection (Version v2.4). Zenodo. http://doi.org/10.5281/zenodo.2678042
#' @param cds the CellDataSet upon which to perform Scrublet
#' @param python_home The python home directory where doubletdetection is installed
#' @param boost_rate (float, optional): Proportion of cell population size to 
#' produce as synthetic doublets.
#' @param n_components (int, optional): Number of principal components used for clustering.
#' @param n_top_var_genes (int, optional): Number of highest variance genes to use; other genes 
#' discarded. Will use all genes when zero. replace (bool, optional): If False, a cell will be 
#' selected as a synthetic doublet's parent no more than once.
#' @param use_phenograph (bool, optional): Set to False to disable PhenoGraph clustering
#' in exchange for louvain clustering implemented in scanpy. Defaults to True.
#' @param n_iters (int, optional): Number of fit operations from which to collect p-values. Defualt value is 25.
#' @param verbose (bool, optional): Set to False to silence all normal operation informational messages. Defaults to True.
#' @param standard_scaling (bool, optional): Set to True to enable standard scaling of normalized count matrix prior to PCA. 
#' Recommended when not using Phenograph. Defaults to False.
#' @param p_thresh (float, optional): hypergeometric test p-value threshold that determines per iteration doublet calls
#' @param voter_thresh (float, optional): fraction of iterations a cell must be called a doublet
#' @return The input CellDataSet with an additional column added to pData with the doublet classification
#' @importFrom reticulate use_python 
#' @importFrom reticulate source_python
#' @references Gayoso, Adam, & Shor, Jonathan. (2018, July 17). DoubletDetection (Version v2.4). Zenodo. http://doi.org/10.5281/zenodo.2678042
#' @export
#' 
doubletdetection<-function(cds, python_home = system("which python", intern = TRUE), module_file=paste(system.file(package = "m3addon"),"doubletdetection.py", sep = "/"),
                           boost_rate=0.25, n_components=30, n_top_var_genes=10000, 
                           use_phenograph=TRUE, n_iters=25,  verbose=TRUE, standard_scaling=FALSE, p_thresh=1e-7, voter_thresh=0.9)
{
  #reticulate::use_python(python_home)
  source_python(module_file)
  X<-t(as.matrix(exprs(cds)))
  doubletdetection_args<-c(list(X, boost_rate=boost_rate, n_components=as.integer(n_components), n_top_var_genes=as.integer(n_top_var_genes), 
                                use_phenograph=use_phenograph, n_iters=as.integer(n_iters), verbose=verbose, 
                                standard_scaling=standard_scaling, p_thresh=p_thresh, voter_thresh=voter_thresh))
  res <- do.call(doubletdetection_py, doubletdetection_args)
  pData(cds)$doubletdetection_call<-res
}
