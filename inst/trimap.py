def trimap_fromR(data, n_dims, n_inliers, n_outliers, n_random, distance, lr, n_iters, knn_tuple, apply_pca, opt_method, verbose, weight_adj, return_seq):
  import trimap
  knn_tuple=None
  embedding = trimap.TRIMAP(int(n_dims), int(n_inliers), int(n_outliers), int(n_random), str(distance), float(lr), int(n_iters), knn_tuple, bool(apply_pca), str(opt_method), bool(verbose), float(weight_adj), bool(return_seq)).fit_transform(data)
  return(embedding)
        
