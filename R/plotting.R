

#' Plot geneset scores
#'
#' @description Geneset scores are a score calculated for each single cell derived from \
#' more than one gene.  This function plots geneset scores using monocle3's 'plot_genes' function.
#' 
#' When using method 'totals', the sum of the size-factor corrected, log-normalized gene \
#' expression for a give set of genes is performed.  When using method 'corrected', single \
#' cell scores for a give gene set that have been "corrected" using 100X genes with similar \
#' expression levels.
 
#' @param cds Input cell_data_set object.
#' @param marker_set Vector of genes in the gene_metadata DataFrame (fData) that can be found in the column 'fData_col'
#' @param name Name given to the geneset
#' @param fData_col Character string denoting the gene_metadata DataFrame (fData) column that contains the genes in marker_set1.  Default = 'gene_short_name'
#' @return Plot
#' @importFrom Matrix colSums
#' @importFrom Matrix t
#' @export
#' @references Puram, S. V. et al. Single-Cell Transcriptomic Analysis of Primary and 
#' Metastatic Tumor Ecosystems in Head and Neck Cancer. Cell 171, 1611.e1–1611.e24 (2017).

plot_geneset<-function(cds, marker_set, name, fData_col="gene_short_name", method=c("totals", "corrected")){
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be one of 'totals' or 'corrected'")
  method <- match.arg(method)
  if(method=="totals") pData(cds)[[name]]<-estimate_score(cds, M2[[name]])
  if(method=="corrected") pData(cds)[[name]]<-estimate_corrected_score(cds, M2[[name]])
  nc<-nchar(name)
  if(nc>50){fontsize<-10}else{fontsize=14}
  switch(method, totals={loca="UW"}, 
         corrected={loca="Broad"})
  plot_cells(cds, color_cells_by = name, label_cell_groups = F, cell_size = 0.5)+ 
    theme(legend.position="top", legend.title = element_blank())+
    ggtitle(paste0(name, ": ", loca))+ 
    theme(plot.title = element_text(size = fontsize, face = "bold"), legend.text = element_text(size=9, angle = 90, vjust=0.5, hjust=0.3))+
    scale_color_gradientn(colors=c( "darkblue","skyblue", "white", "red", "darkred"))
}

monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}




#' Plots the cells along with their trajectories.
#'
#' @param cds cell_data_set for the experiment
#' @param x the column of reducedDims(cds) to plot on the horizontal axis
#' @param y the column of reducedDims(cds) to plot on the vertical axis
#' @param cell_size The size of the point for each cell
#' @param reduction_method The lower dimensional space in which to plot cells.
#'   Must be one of "UMAP", "tSNE", "PCA", "spRing", and "LSI".
#' @param cluster_reduction_method The dimensional space from which to plot clusters.
#'   Must be one of "UMAP", "tSNE", "PCA", "LSI".
#' @param color_cells_by What to use for coloring the cells. Must be either the
#'   name of a column of colData(cds), or one of "clusters", "partitions", or
#'   "pseudotime".
#' @param group_cells_by How to group cells when labeling them. Must be either
#'   the name of a column of colData(cds), or one of "clusters" or "partitions".
#'   If a column in colData(cds), must be a categorical variable.
#' @param genes Facet the plot, showing the expression of each gene in a facet
#'   panel. Must be either a list of gene ids (or short names), or a dataframe
#'   with two columns that groups the genes into modules that will be
#'   aggregated prior to plotting. If the latter, the first column must be gene
#'   ids, and the second must the group for each gene.
#' @param show_trajectory_graph Whether to render the principal graph for the
#'   trajectory. Requires that learn_graph() has been called on cds.
#' @param trajectory_graph_color The color to be used for plotting the
#'   trajectory graph.
#' @param trajectory_graph_segment_size The size of the line segments used for
#'   plotting the trajectory graph.
#' @param norm_method How to normalize gene expression scores prior to plotting
#'   them. Must be one of "log" or "size_only".
#' @param label_cell_groups Whether to label cells in each group (as specified
#'   by group_cells_by) according to the most frequently occurring label(s) (as
#'   specified by color_cells_by) in the group. If false, plot_cells() simply
#'   adds a traditional color legend.
#' @param label_groups_by_cluster Instead of labeling each cluster of cells,
#'   place each label once, at the centroid of all cells carrying that label.
#' @param group_label_size Font size to be used for cell group labels.
#' @param labels_per_group How many labels to plot for each group of cells.
#'   Defaults to 1, which plots only the most frequent label per group.
#' @param label_branch_points Whether to plot a label for each branch point in
#'   the principal graph.
#' @param label_roots Whether to plot a label for each root in the principal
#'   graph.
#' @param label_leaves Whether to plot a label for each leaf node in the
#'   principal graph.
#' @param graph_label_size How large to make the branch, root, and leaf labels.
#' @param alpha Alpha for the cells. Useful for reducing overplotting.
#' @param min_expr Minimum expression threshold for plotting genes
#' @param rasterize Whether to plot cells as a rastered bitmap. Requires the
#'   ggrastr package.
#' @importFrom ggplot2 ggplot
#'
#' @return a ggplot2 plot object
#' @export
#' @examples
#' \dontrun{
#' lung <- load_A549()
#' plot_cells(lung)
#' plot_cells(lung, color_cells_by="log_dose")
#' plot_cells(lung, markers="GDF15")
#' }
#' @references this function differs from that found in monocle3 as it allows for use of spRing \
#' dimensionality reduction.
plot_cells <- function(cds,
                       x=1,
                       y=2,
                       reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "spRing"),
                       color_cells_by="cluster",
                       cluster_reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "spRing"),
                       group_cells_by=c("cluster", "partition"),
                       genes=NULL,
                       show_trajectory_graph=TRUE,
                       trajectory_graph_color="grey28",
                       trajectory_graph_segment_size=0.75,
                       norm_method = c("log", "size_only"),
                       label_cell_groups = TRUE,
                       label_groups_by_cluster=TRUE,
                       group_label_size=2,
                       labels_per_group=1,
                       label_branch_points=TRUE,
                       label_roots=TRUE,
                       label_leaves=TRUE,
                       graph_label_size=2,
                       cell_size=0.35,
                       alpha = 1,
                       min_expr=0.1,
                       rasterize=FALSE) {
  reduction_method <- match.arg(reduction_method)
  cluster_reduction_method <- match.arg(cluster_reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimensions with",
                                      "reduction_method =", reduction_method,
                                      "before attempting to plot."))
  assertthat::assert_that(!is.null(clusters(cds, cluster_reduction_method)),
                          msg = paste("No dimensionality reduction for",
                                      cluster_reduction_method, "calculated.",
                                      "Please run cluster_cells with",
                                      "reduction_method =", reduction_method,
                                      "before attempting to plot."))
  low_dim_coords <- reducedDims(cds)[[reduction_method]]
  assertthat::assert_that(ncol(low_dim_coords) >=max(x,y),
                          msg = paste("x and/or y is too large. x and y must",
                                      "be dimensions in reduced dimension",
                                      "space."))
  if(!is.null(color_cells_by)) {
    assertthat::assert_that(color_cells_by %in% c("cluster", "partition",
                                                  "pseudotime") |
                              color_cells_by %in% names(colData(cds)),
                            msg = paste("color_cells_by must one of",
                                        "'cluster', 'partition', 'pseudotime,",
                                        "or a column in the colData table."))
    
    if(color_cells_by == "pseudotime") {
      tryCatch({pseudotime(cds, reduction_method = reduction_method)},
               error = function(x) {
                 stop(paste("No pseudotime for", reduction_method, "calculated. Please",
                            "run order_cells with reduction_method =", reduction_method,
                            "before attempting to color by pseudotime."))})
      
    }
  }
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers),
                          msg = paste("Either color_cells_by or markers must",
                                      "be NULL, cannot color by both!"))
  
  norm_method = match.arg(norm_method)
  group_cells_by=match.arg(group_cells_by)
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes),
                          msg = paste("Either color_cells_by or genes must be",
                                      "NULL, cannot color by both!"))
  
  if (show_trajectory_graph &&
      is.null(principal_graph(cds)[[reduction_method]])) {
    message("No trajectory to plot. Has learn_graph() been called yet?")
    show_trajectory_graph = FALSE
  }
  
  gene_short_name <- NA
  sample_name <- NA
  #sample_state <- colData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA
  if (rasterize){
    plotting_func <- ggrastr::geom_point_rast
  }else{
    plotting_func <- ggplot2::geom_point
  }
  
  S_matrix <- reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[,c(x,y)])
  
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  
  data_df <- as.data.frame(cbind(data_df, colData(cds)))
  if (group_cells_by == "cluster"){
    data_df$cell_group <-
      tryCatch({clusters(cds,
                         reduction_method = cluster_reduction_method)[
                           data_df$sample_name]},
               error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    data_df$cell_group <-
      tryCatch({partitions(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]},
               error = function(e) {NULL})
  } else{
    stop("Error: unrecognized way of grouping cells.")
  }
  
  if (color_cells_by == "cluster"){
    data_df$cell_color <-
      tryCatch({clusters(cds,
                         reduction_method = cluster_reduction_method)[
                           data_df$sample_name]},
               error = function(e) {NULL})
  } else if (color_cells_by == "partition") {
    data_df$cell_color <-
      tryCatch({partitions(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]},
               error = function(e) {NULL})
  } else if (color_cells_by == "pseudotime") {
    data_df$cell_color <-
      tryCatch({pseudotime(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]}, error = function(e) {NULL})
  } else{
    data_df$cell_color <- colData(cds)[data_df$sample_name,color_cells_by]
  }
  
  ## Graph info
  if (show_trajectory_graph) {
    
    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
      as.data.frame() %>%
      dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
      dplyr::mutate(sample_name = rownames(.),
                    sample_state = rownames(.))
    
    dp_mst <- cds@principal_graph[[reduction_method]]
    
    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      dplyr::select_(source = "from", target = "to") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(
                           source="sample_name",
                           source_prin_graph_dim_1="prin_graph_dim_1",
                           source_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "source") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(
                           target="sample_name",
                           target_prin_graph_dim_1="prin_graph_dim_1",
                           target_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "target")
  }
  
  ## Marker genes
  markers_exprs <- NULL
  expression_legend_label <- NULL
  if (!is.null(genes)) {
    if (!is.null(dim(genes)) && dim(genes) >= 2){
      markers = unlist(genes[,1], use.names=FALSE)
    } else {
      markers = genes
    }
    markers_rowData <- as.data.frame(subset(rowData(cds),
                                            gene_short_name %in% markers |
                                              row.names(rowData(cds)) %in%
                                              markers))
    if (nrow(markers_rowData) == 0) {
      stop("None of the provided genes were found in the cds")
    }
    if (nrow(markers_rowData) >= 1) {
      cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData), ,drop=FALSE]
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))
      
      if (!is.null(dim(genes)) && dim(genes) >= 2){
        genes = as.data.frame(genes)
        row.names(genes) = genes[,1]
        genes = genes[row.names(cds_exprs),]
        
        agg_mat = as.matrix(my.aggregate.Matrix(cds_exprs, as.factor(genes[,2]), fun="sum"))
        
        agg_mat = t(scale(t(log10(agg_mat + 1))))
        agg_mat[agg_mat < -2] = -2
        agg_mat[agg_mat > 2] = 2
        markers_exprs = agg_mat
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        if (is.factor(genes[,2]))
          markers_exprs$feature_id = factor(markers_exprs$feature_id,
                                            levels=levels(genes[,2]))
        
        markers_exprs$feature_label <- markers_exprs$feature_id
        norm_method = "size_only"
        expression_legend_label = "Expression score"
      } else {
        cds_exprs@x = round(cds_exprs@x)
        markers_exprs = matrix(cds_exprs, nrow=nrow(markers_rowData))
        colnames(markers_exprs) = colnames(SingleCellExperiment::counts(cds))
        row.names(markers_exprs) = row.names(markers_rowData)
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        markers_exprs <- merge(markers_exprs, markers_rowData,
                               by.x = "feature_id", by.y="row.names")
        markers_exprs$feature_label <-
          as.character(markers_exprs$gene_short_name)
        
        markers_exprs$feature_label <- ifelse(is.na(markers_exprs$feature_label) | !as.character(markers_exprs$feature_label) %in% markers,
                                              as.character(markers_exprs$feature_id),
                                              as.character(markers_exprs$feature_label))
        
        markers_exprs$feature_label <- factor(markers_exprs$feature_label,
                                              levels = markers)
        if (norm_method == "size_only")
          expression_legend_label = "Expression"
        else
          expression_legend_label = "log10(Expression)"
      }
    }
  }
  
  if (label_cell_groups && is.null(color_cells_by) == FALSE){
    if (is.null(data_df$cell_color)){
      if (is.null(genes)){
        message(paste(color_cells_by, "not found in colData(cds), cells will",
                      "not be colored"))
      }
      text_df = NULL
      label_cell_groups = FALSE
    }else{
      if(is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {
        
        if (label_groups_by_cluster && is.null(data_df$cell_group) == FALSE){
          text_df = data_df %>%
            dplyr::group_by(cell_group) %>%
            dplyr::mutate(cells_in_cluster= dplyr::n()) %>%
            dplyr::group_by(cell_color, add=TRUE) %>%
            dplyr::mutate(per=dplyr::n()/cells_in_cluster)
          median_coord_df = text_df %>%
            dplyr::summarize(fraction_of_group = dplyr::n(),
                             text_x = stats::median(x = data_dim_1),
                             text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df,
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_group) %>%
            dplyr::top_n(labels_per_group, per)
        } else {
          text_df = data_df %>% dplyr::group_by(cell_color) %>%
            dplyr::mutate(per=1)
          median_coord_df = text_df %>%
            dplyr::summarize(fraction_of_group = dplyr::n(),
                             text_x = stats::median(x = data_dim_1),
                             text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df,
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_color) %>%
            dplyr::top_n(labels_per_group, per)
        }
        
        text_df$label = as.character(text_df %>% dplyr::pull(cell_color))
        # I feel like there's probably a good reason for the bit below, but I
        # hate it and I'm killing it for now.
        # text_df$label <- paste0(1:nrow(text_df))
        # text_df$process_label <- paste0(1:nrow(text_df), '_',
        # as.character(as.matrix(text_df[, 1])))
        # process_label <- text_df$process_label
        # names(process_label) <- as.character(as.matrix(text_df[, 1]))
        # data_df[, group_by] <-
        #  process_label[as.character(data_df[, group_by])]
        # text_df$label = process_label
      } else {
        message(paste("Cells aren't colored in a way that allows them to",
                      "be grouped."))
        text_df = NULL
        label_cell_groups = FALSE
      }
    }
  }
  
  if (!is.null(markers_exprs) && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name",
                     by.y="cell_id")
    data_df$value <- with(data_df, ifelse(value >= min_expr, value, NA))
    na_sub <- data_df[is.na(data_df$value),]
    if(norm_method == "size_only"){
      g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
        plotting_func(aes(data_dim_1, data_dim_2), size=I(cell_size),
                      stroke = I(cell_size / 2), color = "grey80",
                      data = na_sub) +
        plotting_func(aes(color=value), size=I(cell_size),
                      stroke = I(cell_size / 2), na.rm = TRUE) +
        viridis::scale_color_viridis(option = "viridis",
                                     name = expression_legend_label,
                                     na.value = "grey80", end = 0.8) +
        guides(alpha = FALSE) + facet_wrap(~feature_label)
    } else {
      g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
        plotting_func(aes(data_dim_1, data_dim_2), size=I(cell_size),
                      stroke = I(cell_size / 2), color = "grey80",
                      data = na_sub) +
        plotting_func(aes(color=log10(value+min_expr)),
                      size=I(cell_size), stroke = I(cell_size / 2),
                      na.rm = TRUE) +
        viridis::scale_color_viridis(option = "viridis",
                                     name = expression_legend_label,
                                     na.value = "grey80", end = 0.8) +
        guides(alpha = FALSE) + facet_wrap(~feature_label)
    }
  } else {
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))
    
    # We don't want to force users to call order_cells before even being able
    # to look at the trajectory, so check whether it's null and if so, just
    # don't color the cells
    if(color_cells_by %in% c("cluster", "partition")){
      if (is.null(data_df$cell_color)){
        g <- g + geom_point(color=I("gray"), size=I(cell_size), na.rm = TRUE,
                            alpha = I(alpha))
        message(paste("cluster_cells() has not been called yet, can't",
                      "color cells by cluster"))
      } else{
        g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                            na.rm = TRUE, alpha = alpha)
      }
      g <- g + guides(color = guide_legend(title = color_cells_by,
                                           override.aes = list(size = 4)))
    } else if (class(data_df$cell_color) == "numeric"){
      g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                          na.rm = TRUE, alpha = alpha)
      g <- g + viridis::scale_color_viridis(name = color_cells_by, option="C")
    } else {
      g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                          na.rm = TRUE, alpha = alpha)
      g <- g + guides(color = guide_legend(title = color_cells_by,
                                           override.aes = list(size = 4)))
    }
    
  }
  if (show_trajectory_graph){
    g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1",
                                     y="source_prin_graph_dim_2",
                                     xend="target_prin_graph_dim_1",
                                     yend="target_prin_graph_dim_2"),
                          size=trajectory_graph_segment_size,
                          color=I(trajectory_graph_color),
                          linetype="solid",
                          na.rm=TRUE,
                          data=edge_df)
    
    
    if (label_branch_points){
      mst_branch_nodes <- branch_nodes(cds)
      branch_point_df <- ica_space_df %>%
        dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
        dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
      
      g <- g +
        geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                   shape = 21, stroke=I(trajectory_graph_segment_size),
                   color="white",
                   fill="black",
                   size=I(graph_label_size * 1.5),
                   na.rm=TRUE, branch_point_df) +
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                             label="branch_point_idx"),
                  size=I(graph_label_size), color="white", na.rm=TRUE,
                  branch_point_df)
    }
    
    if (label_leaves){
      mst_leaf_nodes <- leaf_nodes(cds)
      leaf_df <- ica_space_df %>%
        dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>%
        dplyr::mutate(leaf_idx = seq_len(dplyr::n()))
      
      g <- g +
        geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                   shape = 21, stroke=I(trajectory_graph_segment_size),
                   color="black",
                   fill="lightgray",
                   size=I(graph_label_size * 1.5),
                   na.rm=TRUE,
                   leaf_df) +
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                             label="leaf_idx"),
                  size=I(graph_label_size), color="black", na.rm=TRUE, leaf_df)
    }
    
    if (label_roots){
      mst_root_nodes <- root_nodes(cds)
      root_df <- ica_space_df %>%
        dplyr::slice(match(names(mst_root_nodes), sample_name)) %>%
        dplyr::mutate(root_idx = seq_len(dplyr::n()))
      
      g <- g +
        geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                   shape = 21, stroke=I(trajectory_graph_segment_size),
                   color="black",
                   fill="white",
                   size=I(graph_label_size * 1.5),
                   na.rm=TRUE,
                   root_df) +
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                             label="root_idx"),
                  size=I(graph_label_size), color="black", na.rm=TRUE, root_df)
    }
  }
  
  if(label_cell_groups) {
    g <- g + ggrepel::geom_text_repel(data = text_df,
                                      mapping = aes_string(x = "text_x",
                                                           y = "text_y",
                                                           label = "label"),
                                      size=I(group_label_size))
    # If we're coloring by gene expression, don't hide the legend
    if (is.null(markers_exprs))
      g <- g + theme(legend.position="none")
  }
  
  g <- g +
    #scale_color_brewer(palette="Set1") +
    monocle_theme_opts() +
    xlab(paste(reduction_method, x)) +
    ylab(paste(reduction_method, y)) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white'))
  g
}

#' GSEA and plot
#' @description Performs GSEA (using fgsea) and returns a GSEA-style enrichment plot
#' @param pathway = vector of genes
#' @param stats = vectir of ranked gene stats (usually fold change or SNR) with names \
#' that contain the gene name
#' @param rug = whether to make a rug plot
#' @param dot_enhance character string denoting a color that enhances the dot appearance \
#'  with another color
#' @import ggplot2
#' @importFrom fgsea calcGseaStat
#' @param all_the_rest_of_them Should be self explanatory
#' @return Performs GSEA of "pathway" genes on stats'
#' @references fgsea package
#' @export


enrichmentPlot<-function (pathway, stats, 
                          gseaParam = 1, 
                          segment=F, rug=T, 
                          rug_color="black", segment_color="black", 
                          dot_color="green", dot_enhance=NULL, 
                          dot_enhance_size=2, dot_shape=21, 
                          dot_enhance_alpha=0.7, dot_size=1,
                          return_data=FALSE,
                          print_plot=FALSE,
                          return_plot=TRUE) 
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE, returnLeadingEdge = TRUE)
  #bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(pathway)
  ys <- as.vector(rbind(tops))
  le <- c(rep(1, length(xs)))
  le[xs %in% gseaRes$le]<-5
  le_bool<-le==5
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0), le=c(0,le,0))
  diff <- (max(ys) - min(ys))/6
  df_out<-data.frame(Rank = c(xs), ES = c(ys), LE=le_bool, row.names = names(tops))
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) + geom_point(color = dot_color, 
                                                      size = dot_size) + geom_hline(yintercept = max(ys), colour = "black", 
                                                                                    linetype = "dashed") + geom_hline(yintercept = min(ys), 
                                                                                                                      colour = "black", linetype = "dashed") + geom_hline(yintercept = 0, 
                                                                                                                                                                          colour = "black") + geom_line(color = segment_color) + theme_bw()
  # if(rug) {g<-g+geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, 
  #                                                                          y = min(ys)-diff/2, xend = x, yend = min(ys)-diff), size = 0.2, color=rug_color)}  
  if(!is.null(dot_enhance)) {g<-g+geom_point(color = dot_enhance,
                                             size = dot_enhance_size, shape=dot_shape, alpha=dot_enhance_alpha)}
  
  if(rug){ g<-g+ geom_point(data = toPlot, aes(x = x,
                                               y = min(ys)-diff, size = le), colour = rug_color,shape = 124)}
  g<-g+theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
    labs(x = "rank", y = "enrichment score")+ theme(legend.position="none")
  if(print_plot) print(g)
  if(return_data) return(list(plot=g, gseaRes=gseaRes, df_out=df_out))
  if(return_plot) return(g)
}

#' @keywords internal
compileStats<-function(gsa, gs=NULL){
  #this function collapses adjusted stats from a gsa (piano package) generated using gsea or fgsea
  stats<-data.table::data.table(up=as.vector(gsa$pAdjDistinctDirUp), down=as.vector(gsa$pAdjDistinctDirDn))
  if(is.null(gs)){
    gs<-names(gsa$gsc)
  }
  stats[is.na(stats)]<-1
  stats<-1-stats
  dir<-colnames(stats)[max.col(stats, ties.method = "first")]
  stats<-abs(stats-1)
  stats$dir<-dir
  stats$name<-names(gsa$gsc)
  data.table::setkey(data.table::as.data.table(stats), name)
  return( stats[gs])
}

#' @export
returnFDR<-function(gsa, gs=NULL){
  #this function returns a stat line describing the FDR and direction of a gsa generated using gsea or fsea
  dat<-compileStats(gsa=gsa, gs=gs)
  paste0("FDR = ", round(dat[[dat$dir]], 6), " in the ",dat$dir,  " direction")
}

#' @keywords internal
monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}


