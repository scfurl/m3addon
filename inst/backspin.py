
def backSPIN (data, numLevels, first_run_iters, first_run_step, runs_iters, runs_step, \
split_limit_g, split_limit_c, stop_const, low_thrs, verbose):
  #print(locals())
  from backspinpy import SPIN, backSPIN, fit_CV, feature_selection, CEF_obj
  
  bs = backSPIN(data=data, numLevels=int(numLevels), first_run_iters=int(first_run_iters), \
  first_run_step=first_run_step, runs_iters=int(runs_iters), runs_step=runs_step, \
  split_limit_g=int(split_limit_g), split_limit_c=int(split_limit_c), stop_const = stop_const, \
  low_thrs=low_thrs, verbose=verbose)
  
  return(bs.genes_order, bs.cells_order, bs.genes_gr_level, bs.cells_gr_level, \
  bs.cells_gr_level_sc, bs.genes_bor_level, bs.cells_bor_level)

'''Run the backSPIN algorithm
Parameters
----------
data: 2-D array
the data matrix, rows should be genes and columns single cells/samples
numLevels: int
the number of splits that will be tried
first_run_iters: float
the iterations of the preparatory SPIN
first_run_step: float
the step parameter passed to _generate_widlist for the preparatory SPIN
runs_iters: int
the iterations parameter passed to the _divide_to_2and_resort.
influences all the SPIN iterations except the first
runs_step: float
the step parameter passed to the _divide_to_2and_resort.
influences all the SPIN iterations except the first
wid: float
the wid of every iteration of the splitting and resorting
split_limit_g: int
If the number of specific genes in a subgroup is smaller than this number
splitting of that subgrup is not allowed
split_limit_c: int
If the number cells in a subgroup is smaller than this number splitting of
that subgrup is not allowed
stop_const: float
minimum score that a breaking point has to reach to be suitable for splitting
low_thrs: float
genes with average lower than this threshold are assigned to either of the
splitting group reling on genes that are higly correlated with them

Returns
-------
results: Result object
The results object contain the following attributes
genes_order: 1-D array
  indexes (a permutation) sorting the genes
cells_order: 1-D array
  indexes (a permutation) sorting the cells
genes_gr_level: 2-D array
  for each depth level contains the cluster indexes for each gene
cells_gr_level:
  for each depth level contains the cluster indexes for each cell
cells_gr_level_sc:
  score of the splitting
genes_bor_level:
  the border index between gene clusters
cells_bor_level:
  the border index between cell clusters

Notes
-----
Typical usage
'''
