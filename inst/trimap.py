def trimap_fromR(data, n_dims, n_inliers, n_outliers, n_random, distance, lr, n_iters, knn_tuple, apply_pca, opt_method, verbose, weight_adj, return_seq):
  import trimap
  knn_tuple=None
  embedding = trimap.TRIMAP(n_dims = int(n_dims), n_inliers = int(n_inliers), n_outliers = int(n_outliers), n_random = int(n_random), distance = str(distance), lr = float(lr), n_iters = int(n_iters), apply_pca = bool(apply_pca),opt_method = str(opt_method), verbose = bool(verbose), weight_adj = float(weight_adj), return_seq = bool(return_seq)).fit_transform(data)
  return(embedding)
        
