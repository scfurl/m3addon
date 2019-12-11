

#' Backspin
#' @description The BackSPIN biclustering algorithm was developed by Amit Zeisel and is described in Zeisel et al. 
#' Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq Science 2015 (PMID: 25700174, 
#' doi: 10.1126/science.aaa1934). Please cite this paper if you use the BackSPIN algorithm in your work.
#' Original MATLAB implementation by Amit Zeisel. This repo contains a standalone command-line version of 
#' BackSPIN, implemented in Python by Gioele La Manno.
#' @param data Rows should be genes and columns single cells/samples.
#' @param levels The number of nested splits that will be tried by the algorithm
#' @param first_run_iters Number of the iterations used in the preparatory SPIN. Default 10
#' @param first_run_step the step parameter passed to _generate_widlist for the preparatory SPIN. Default 0.05
#' @param runs_iters the iterations parameter passed to the _divide_to_2and_resort. influences all the SPIN iterations except the first
#' @param runs_step: the step parameter passed to the _divide_to_2and_resort. influences all the SPIN iterations except the first
#' @param split_limit_g: If the number of specific genes in a subgroup is smaller than this number splitting of that subgrup is not allowed
#' @param split_limit_c: If the number cells in a subgroup is smaller than this number splitting of that subgrup is not allowed
#' @param stop_const: minimum score that a breaking point has to reach to be suitable for splitting
#' @param low_thrs: genes with average lower than this threshold are assigned to either of the splitting group reling on genes that are higly correlated with them
#' @param verbose Print the extra details of what is happening
#' @export
#' 
#' @importFrom reticulate source_python py_available
#' @references Zeisel et. al. Science. 2015 Mar 6;347(6226):1138-42. doi: 10.1126/science.aaa1934. Epub 2015 Feb 19.

backspin<-function(data, numLevels=2, first_run_iters=10.0, first_run_step=0.05, runs_iters=8 ,runs_step=0.25,
                   split_limit_g=2, split_limit_c=2, stop_const = 1.15, low_thrs=0.2, verbose=T){
  if(!py_available("backspin")) stop("python module backspin does not seem to be installed; - try running 'py_config()'")
  source_python(paste(system.file(package = "m3addon"), 
                                  "backspin.py", sep = "/"))
  
  backspin_args<-c(list(data=data, numLevels=numLevels, first_run_iters=first_run_iters, first_run_step=first_run_step,
                        runs_iters=runs_iters, runs_step=runs_step, split_limit_g=split_limit_g, split_limit_c=split_limit_c,
                        stop_const = stop_const, low_thrs=low_thrs, verbose=verbose))
  bs = do.call(backSPIN, backspin_args)
  names(bs) = unlist(strsplit("genes_order, cells_order, genes_gr_level, cells_gr_level, cells_gr_level_sc, genes_bor_level, cells_bor_level", "\\,"))
  return(bs)
}
