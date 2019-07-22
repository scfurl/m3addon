#' @export

# backspin<-function(data, numLevels=2, first_run_iters=10.0, first_run_step=0.05, runs_iters=8 ,runs_step=0.25,
#                    split_limit_g=2, split_limit_c=2, stop_const = 1.15, low_thrs=0.2, verbose=F){
#   reticulate::source_python(paste(system.file(package = "m3addon"), 
#                                   "backspin.py", sep = "/"))
#   backspin_args<-c(list(data=data, numLevels=numLevels, first_run_iters=first_run_iters, first_run_step=first_run_step,
#   runs_iters=runs_iters, runs_step=runs_step, split_limit_g=split_limit_g, split_limit_c=split_limit_c,
#   stop_const = stop_const, low_thrs=low_thrs, verbose=verbose))
#   do.call(backSPIN, backspin_args)
# }

backspin<-function(data, numLevels=2, first_run_iters=10.0, first_run_step=0.05, runs_iters=8 ,runs_step=0.25,
                   split_limit_g=2, split_limit_c=2, stop_const = 1.15, low_thrs=0.2, verbose=T){
  reticulate::source_python(paste(system.file(package = "m3addon"), 
                                  "backspin.py", sep = "/"))
  backspin_args<-c(list(data=data, numLevels=numLevels, first_run_iters=first_run_iters, first_run_step=first_run_step,
                        runs_iters=runs_iters, runs_step=runs_step, split_limit_g=split_limit_g, split_limit_c=split_limit_c,
                        stop_const = stop_const, low_thrs=low_thrs, verbose=verbose))
  bs = do.call(backSPIN, backspin_args)
  names(bs) = unlist(strsplit("genes_order, cells_order, genes_gr_level, cells_gr_level, cells_gr_level_sc, genes_bor_level, cells_bor_level", "\\,"))
  return(bs)
}