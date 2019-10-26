#' @export
plot_gene_dispersion<-function(cds){
  if(is.null(cds@int_metadata$dispersion$disp_func)){
    prd<-cds@int_metadata$dispersion
    g<-ggplot2::ggplot(prd, ggplot2::aes(x = log_mean, y = fit)) 
    if("use_for_ordering" %in% colnames(cds@int_metadata$dispersion)){
      g <- g + ggplot2::geom_point(data=prd, ggplot2::aes(x=log_mean, y=log_dispersion, color=use_for_ordering), alpha=0.4)
    }else{
      g <- g + ggplot2::geom_point(data=prd, ggplot2::aes( x=log_mean, y=log_dispersion, color="grey"), alpha=0.4)
    }
    g<-g+
      ggplot2::theme_bw() +
      ggplot2::geom_line() # + 
      #ggplot2::geom_smooth(data=prd, ggplot2::aes(ymin = lci, ymax = uci), stat = "identity")
    return(g)
  }else{
    prd<-cds@int_metadata$dispersion$disp_table
    prd$fit<-log(cds@int_metadata$dispersion$disp_func(prd$mu))
    prd$mu<-log(prd$mu)
    prd$disp<-log(prd$disp)
    colnames(prd)<-c("gene_id", "log_mean", "log_dispersion", "fit")
    g<-ggplot2::ggplot(prd, ggplot2::aes(x = log_mean, y = fit)) 
    if("use_for_ordering" %in% names(cds@int_metadata$dispersion)){
      prd$use_for_ordering = cds@int_metadata$dispersion$use_for_ordering
      g <- g + ggplot2::geom_point(data=prd, ggplot2::aes(x=log_mean, y=log_dispersion, color=use_for_ordering, alpha=0.4))
    }else{
      g <- g + ggplot2::geom_point(data=prd, ggplot2::aes( x=log_mean, y=log_dispersion, color="grey", alpha=0.4))
    }
    g<-g+
      ggplot2::theme_bw() +
      ggplot2::geom_line(data=prd, ggplot2::aes( x=log_mean, y=fit)) +
      ggplot2::geom_smooth(data=prd, formula = fit ~ log_mean, stat = "identity")
    return(g)
  }
}

#' Select genes in a cell_data_set for dimensionality reduction
#'
#' @description Monocle3 aims to learn how cells transition through a
#' biological program of gene expression changes in an experiment. Each cell
#' can be viewed as a point in a high-dimensional space, where each dimension
#' describes the expression of a different gene. Identifying the program of
#' gene expression changes is equivalent to learning a \emph{trajectory} that
#' the cells follow through this space. However, the more dimensions there are
#' in the analysis, the harder the trajectory is to learn. Fortunately, many
#' genes typically co-vary with one another, and so the dimensionality of the
#' data can be reduced with a wide variety of different algorithms. Monocle3
#' provides two different algorithms for dimensionality reduction via
#' \code{reduce_dimensions} (UMAP and tSNE). The function
#' \code{select_features} is an optional step in the trajectory building
#' process before \code{preprocess_cds}.  After calculating dispersion for
#' a cell_data_set using the \code{calculate_gene_dispersion} function, the 
#' \code{select_genes} function allows the user to identify a set of genes
#' that will be used in downstream dimensionality reduction methods.
#'
#'
#' @param cds the cell_data_set upon which to perform this operation.
#' @param fit_min the minimum multiple of the dispersion fit calculation; default = 1
#' @param fit_max the maximum multiple of the dispersion fit calculation; default = Inf
#' @param logmean_ul the maximum multiple of the dispersion fit calculation; default = Inf
#' @param logmean_ll the maximum multiple of the dispersion fit calculation; default = Inf
#' @return an updated cell_data_set object with selected features 
#' @export

select_genes<-function(cds, fit_min=1, fit_max=Inf, logmean_ul=NULL, logmean_ll=NULL){
  if(is.null(cds@int_metadata$dispersion$disp_func)){
    df<-cds@int_metadata$dispersion
    df$ratio<-df$log_dispersion/df$fit
    cds@int_metadata$dispersion$use_for_ordering <- df$ratio > fit_min & df$ratio < fit_max
    if(!is.null(logmean_ll)){
      cds@int_metadata$dispersion$use_for_ordering <- cds@int_metadata$dispersion$use_for_ordering & df$log_mean > logmean_ll
    }
    if(!is.null(logmean_ul)){
      cds@int_metadata$dispersion$use_for_ordering <- cds@int_metadata$dispersion$use_for_ordering & df$log_mean < logmean_ul
    }
    return(cds)
  }else{
    df<-cds@int_metadata$dispersion$disp_table
    df$fit<-cds@int_metadata$dispersion$disp_func(df$mu)
    df$ratio<-df$disp/df$fit
    df$log_disp=log(df$disp)
    df$log_mean<-log(df$mu)
    cds@int_metadata$dispersion$use_for_ordering <- df$ratio > fit_min & df$ratio < fit_max
    if(!is.null(logmean_ll)){
      cds@int_metadata$dispersion$use_for_ordering <- cds@int_metadata$dispersion$use_for_ordering & df$log_mean > logmean_ll
    }
    if(!is.null(logmean_ul)){
      cds@int_metadata$dispersion$use_for_ordering <- cds@int_metadata$dispersion$use_for_ordering & df$log_mean < logmean_ul
    }
    return(cds)
  }
}


#' @export
get_ordering_genes<-function(cds, gene_column="id"){
  if(is.null(cds@int_metadata$dispersion$disp_func)){
    as.character(cds@int_metadata$dispersion[[gene_column]][cds@int_metadata$dispersion$use_for_ordering])
  }else{
    as.character(cds@int_metadata$dispersion$disp_table[[gene_column]][cds@int_metadata$dispersion$use_for_ordering])
  }
}

#' Calculate dispersion genes in a cell_data_set object
#'
#' @description Monocle3 aims to learn how cells transition through a
#' biological program of gene expression changes in an experiment. Each cell
#' can be viewed as a point in a high-dimensional space, where each dimension
#' describes the expression of a different gene. Identifying the program of
#' gene expression changes is equivalent to learning a \emph{trajectory} that
#' the cells follow through this space. However, the more dimensions there are
#' in the analysis, the harder the trajectory is to learn. Fortunately, many
#' genes typically co-vary with one another, and so the dimensionality of the
#' data can be reduced with a wide variety of different algorithms. Monocle3
#' provides two different algorithms for dimensionality reduction via
#' \code{reduce_dimensions} (UMAP and tSNE). The function
#' \code{calculate_dispersion} is an optional step in the trajectory building
#' process before \code{preprocess_cds}.  After calculating dispersion for
#' a cell_data_set using the \code{calculate_gene_dispersion} function, the 
#' \code{select_genes} function allows the user to identify a set of genes
#' that will be used in downstream dimensionality reduction methods.  These
#' genes and their disperion and mean expression can be plotted using the 
#' \code{plot_gene_dispersion} function.
#'
#'
#' @param cds the cell_data_set upon which to perform this operation.
#' @param q the polynomial degree; default = 3.
#' @param id_tag the name of the feature data column corresponding to 
#' the unique id - typically ENSEMBL id; default = "id".
#' @param symbol_tag the name of the feature data column corresponding to 
#' the gene symbol; default = "gene_short_name".
#' @return an updated cell_data_set object with dispersion and mean expression saved
#' @export
calculate_gene_dispersion<-function(cds, q=3, id_tag="id", symbol_tag="gene_short_name", method, removeOutliers=T){
  if(method=="m2"){
    df<-data.frame(calc_dispersion_m2(cds = cds, min_cells_detected = 1, min_exprs=1))
    fdat<-fData(cds)
    if (!is.list(df)) 
      stop("Parametric dispersion fitting failed, please set a different lowerDetectionLimit")
    disp_table <- subset(df, is.na(mu) == FALSE)
    res <- monocle:::parametricDispersionFit(disp_table, verbose = T)
    fit <- res[[1]]
    coefs <- res[[2]]
    if (removeOutliers) {
        CD <- cooks.distance(fit)
        cooksCutoff <- 4/nrow(disp_table)
        message(paste("Removing", length(CD[CD > cooksCutoff]), 
                      "outliers"))
        outliers <- union(names(CD[CD > cooksCutoff]), setdiff(row.names(disp_table), 
                                                               names(CD)))
        res <- monocle:::parametricDispersionFit(disp_table[row.names(disp_table) %in% 
                                                    outliers == FALSE, ], verbose=T)
        fit <- res[[1]]
        coefs <- res[[2]]
      names(coefs) <- c("asymptDisp", "extraPois")
      ans <- function(q) coefs[1] + coefs[2]/q
      attr(ans, "coefficients") <- coefs
    }
    res <- list(disp_table = disp_table, disp_func = ans)
    cds@int_metadata$dispersion<-res
    return(cds)
  }
  
  if(method=="m2like"){
    rounded <- round(exprs(cds))
    nzGenes <- Matrix::rowSums(rounded > min_exprs)
    nzGenes <- names(nzGenes[nzGenes > min_cells_detected])
    x <- DelayedArray(t(t(rounded[nzGenes, ])/pData(cds[nzGenes, 
                                                        ])$Size_Factor))
    xim <- mean(1/pData(cds[nzGenes, ])$Size_Factor)
    if (class(exprs(cds)) %in% c("dgCMatrix", "dgTMatrix")) {
      f_expression_mean <- as(DelayedMatrixStats::rowMeans2(x), 
                              "sparseVector")
    }
    else {
      f_expression_mean <- DelayedMatrixStats::rowMeans2(x)
    }
    f_expression_var <- DelayedMatrixStats::rowVars(x)
    # m3vars<-m3addon:::rowStdDev(as(x, "sparseMatrix"))
    # head(m3vars[1,]^(2))
    # head(f_expression_var)
    disp_guess_meth_moments <- f_expression_var - xim * f_expression_mean
    disp_guess_meth_moments <- disp_guess_meth_moments/(f_expression_mean^2)
    res <- data.frame(mu = as.vector(f_expression_mean), disp = as.vector(disp_guess_meth_moments))
    res[res$mu == 0]$mu = NA
    res[res$mu == 0]$disp = NA
    res$disp[res$disp < 0] <- 0
    res$mu<-log(res$mu)
    res$disp<-log(res$disp)
    res <- cbind(gene_id = row.names(fData(cds[nzGenes, ])), 
                 res)
    res<-res[is.finite(res$disp),]
    model <- lm(data = res, disp ~ mu + poly(mu, degree=q))
    prd <- data.frame(mu = res$mu, id=res$gene_id)
    err<-suppressWarnings(predict(model, newdata= prd, se.fit = T))
    prd$log_mean<-prd$mu
    prd$lci <- err$fit - 1.96 * err$se.fit
    prd$fit <- err$fit
    prd$uci <- err$fit + 1.96 * err$se.fit
    prd$log_dispersion<-res$disp
    prd$id<-res[["gene_id"]]
    prd$gene_short_name<-fdat[[symbol_tag]][match(prd$id, fdat[[id_tag]])]
    cds@int_metadata$dispersion<-prd
    return(cds)
  }
  
  if(method=="m3addon"){
    m<-Matrix::rowMeans(counts(cds))
    # sd<-sqt(rowVars(as.matrix(counts(cds))))
    sd<-m3addon:::rowStdDev(exprs(cds))[1,]
    fdat<-fData(cds)
    cv<-sd/m*100
    df<-data.frame(log_dispersion=log(cv), log_mean=log(m))
    df[[id_tag]]<-fdat[[id_tag]]
    df<-df[is.finite(df$log_dispersion),]
    model <- lm(data = df, log_dispersion ~ log_mean + poly(log_mean, degree=q))
    prd <- data.frame(log_mean = df$log_mean)
    err<-suppressWarnings(predict(model, newdata= prd, se.fit = T))
    prd$lci <- err$fit - 1.96 * err$se.fit
    prd$fit <- err$fit
    prd$uci <- err$fit + 1.96 * err$se.fit
    prd$log_dispersion<-df$log_dispersion
    prd$id<-df[[id_tag]]
    prd$gene_short_name<-fdat[[symbol_tag]][match(prd$id, fdat[[id_tag]])]
    cds@int_metadata$dispersion<-prd
    return(cds)
  }
}

calc_dispersion_m2<-function (cds, expressionFamily, min_cells_detected=1, min_exprs = 1) 
{
  rounded <- round(exprs(cds))
  nzGenes <- Matrix::rowSums(rounded > min_exprs)
  nzGenes <- names(nzGenes[nzGenes > min_cells_detected])
  x <- DelayedArray(t(t(rounded[nzGenes, ])/pData(cds[nzGenes, 
                                                      ])$Size_Factor))
  xim <- mean(1/pData(cds[nzGenes, ])$Size_Factor)
  if (class(exprs(cds)) %in% c("dgCMatrix", "dgTMatrix")) {
    f_expression_mean <- as(DelayedMatrixStats::rowMeans2(x), 
                            "sparseVector")
  }
  else {
    f_expression_mean <- DelayedMatrixStats::rowMeans2(x)
  }
  f_expression_var <- DelayedMatrixStats::rowVars(x)
  # m3vars<-m3addon:::rowStdDev(as(x, "sparseMatrix"))
  # head(m3vars[1,]^(2))
  # head(f_expression_var)
  disp_guess_meth_moments <- f_expression_var - xim * f_expression_mean
  disp_guess_meth_moments <- disp_guess_meth_moments/(f_expression_mean^2)
  res <- data.frame(mu = as.vector(f_expression_mean), disp = as.vector(disp_guess_meth_moments))
  res[res$mu == 0]$mu = NA
  res[res$mu == 0]$disp = NA
  res$disp[res$disp < 0] <- 0
  res <- cbind(gene_id = row.names(fData(cds[nzGenes, ])), 
               res)
  res
}











    

  
