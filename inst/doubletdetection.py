def doubletdetection_py(X, boost_rate, n_components, n_top_var_genes, use_phenograph, n_iters, normalizer, verbose, standard_scaling, \
                        p_thresh, voter_thresh):
  import doubletdetection
  if use_phenograph:
    phenograph_parameters = {"prune":True}
  clf = doubletdetection.BoostClassifier(boost_rate, int(n_components), int(n_top_var_genes), use_phenograph, \
                        phenograph_parameters, int(n_iters), normalizer, verbose, standard_scaling)
  # raw_counts is a cells by genes count matrix
  labels = clf.fit(X).predict(p_thresh, voter_thresh)
  return(labels)
