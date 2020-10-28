#'  Still working on this function
#'  @export
cluster_PCA<-function(cds, method="louvain",
                      k = 20, dims=NULL,
                      weight = F, 
                      num_iter = 1, 
                      resolution_parameter = NULL, 
                      random_seed = 2020, 
                      verbose = T, 
                      partition_q_value = 0.05)
{
  if(is.null(dims)){
    dims<-1:dim(reducedDims(cds)$PCA)[2]
  }
  if(method=="leiden"){
    cmethod<-monocle3:::leiden_clustering
  }
  else{
    cmethod<-monocle3:::louvain_clustering
  }
  cluster_result<-cmethod(data = as.matrix(reducedDims(cds)$PCA)[,dims], 
                          pd = colData(cds), k = k, weight = weight, 
                          num_iter = num_iter, 
                          resolution_parameter = resolution_parameter, 
                          random_seed = random_seed, verbose = verbose)
  cluster_graph_res<-monocle3:::compute_partitions(cluster_result$g, 
                                                   cluster_result$optim_res, partition_q_value, verbose=t)
  partitions <- as.factor(igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership])
  clusters <- factor(igraph::membership(cluster_result$optim_res))
  cds@clusters[["UMAP"]] <- list(cluster_result = cluster_result, 
                                 partitions = partitions, clusters = clusters)
  cds
}
