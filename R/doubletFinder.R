#' doubletFinder
#'
#' @description Core doublet prediction function of the DoubletFinder package. Generates artifical 
#' doublets from an existing, pre-processed cdsrat object. Real and artificial data are then merged 
#' and pre-processed using parameters utilized for the existing cdsrat object. PC distance matrix is 
#' then computed and used the measure the proportion of artificial nearest neighbors (pANN) for every 
#' real cell. pANN is then thresholded according to the number of expected doublets to generate final 
#' doublet predictions.
#' @usage cds <- doubletFinder(cds, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE)
#' @param cds Input cell_data_set object.
#' @param PCs	Number of statistically-significant principal components (e.g., as estimated from PC 
#' elbow plot)
#' @param ...  arguments passed to 1)  \code{calculate_gene_dispersion} and 2) \code{select_genes}: 
#' note that default for top_n passsed to \code{select_genes} in this context is 2000 features.  See specific function documentation for further information
#' on acceptable arguments to be passed these functions.
#' @param pN	The number of generated artificial doublets, expressed as a proportion of the merged 
#' real-artificial data. Default is set to 0.25, based on observation that DoubletFinder performance is
#'  largely pN-invariant (see McGinnis, Murrow and Gartner 2019, Cell Systems).
#' @param pK The PC neighborhood size used to compute pANN, expressed as a proportion of the merged 
#' real-artificial data. No default is set, as pK should be adjusted for each scRNA-seq dataset. Optimal 
#' pK values can be determined using mean-variance-normalized bimodality coefficient.
#' @param nExp The total number of doublet predictions produced. This value can best be estimated 
#' from cell loading densities into the 10X/Drop-Seq device, and adjusted according to the estimated 
#' proportion of homotypic doublets.
#' @param genes if "all", use all genes; if "recalc", recalculate ordering genes using \code{calculate_gene_dispersion} and 
#' \code{select_genes}, passing arguments to each of these functions using ....; if "same" use ordering genes 
#' specified in cds; or a vector of ordering genes to be used.
#' @param reuse.pANN Metadata column name for previously-generated pANN results. Argument should be set to 
#' FALSE (default) for initial DoubletFinder runs. Enables fast adjusting of doublet predictions for 
#' different nExp.
#' @importFrom fields rdist
#' @references McGinnis, Murrow and Gartner 2019, Cell Systems; https://github.com/chris-mcginnis-ucsf/DoubletFinder
#' @export


doubletFinder_v3 <- function(cds, PCs=1:100, pN = 0.25, pK, nExp, genes=c("same", "all", "recalc"), ...) {
    dots <- list(...)
    og_done<-FALSE
    sg_args <-c("logmean_ul", "logmean_ll", "top_n", "fit_min", "fit_max")
    cd_args<-c("id", "symbol_tag","method", "remove_outliers", "q")
    #gather relevant args
    rel_args<-dots[names(dots) %in% c(sg_args, cd_args)]
    if(!"top_n" %in% names(rel_args)){
      rel_args<-c(list(top_n=2000), rel_args)
    }
    if(length(genes)>1 & all(genes %in% c("same", "all", "recalc"))){
      genes <- "same"
    }
    if(length(genes)==1 & genes %in% c("same", "all", "recalc")){
      if(genes=="all") {
        message("Using all features")
        ord_genes <- rownames(fData(cds))
        og_done <-TRUE
      }
      if(genes=="recalc"){
        message("Recalculating ordering features using the following arguments:\nCalculate Dispersion:\n")
        print(rel_args[names(rel_args) %in% cd_args])
        message("\nSelect Genes:\n")
        print(rel_args[names(rel_args) %in% sg_args])
        rel_args<-c(list(cds=cds), rel_args)
        cds<-do.call(calculate_gene_dispersion, rel_args[names(rel_args) %in% c("cds", cd_args)])
        cds<-do.call(select_genes, rel_args[names(rel_args) %in% c("cds", sg_args)])
        ord_genes<-get_ordering_genes(cds)
        og_done <-TRUE
      }
      if(genes=="same"){
        message("Using existing ordering features")
        ord_genes<-get_ordering_genes(cds)
        og_done <-TRUE
      }
    }
    if(length(genes)>1 & all(genes %in% rownames(exprs(cds)))){
      ord_genes<-genes
      message("Using supplied ordering genes")
      og_done <-TRUE
    }
    if(!og_done & all(genes %in% rownames(exprs(cds)))){
      stop("Genes not found in cds; they must be rownames of exprs(cds)")
    }
    real.cells <- rownames(cds@colData)
    data <- exprs(cds)[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    message(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    
    ## Pre-process cds object
    message("Creating Monocle3 object with doublets...")
    cds_wdoublets <- new_cell_data_set(data_wdoublets, gene_metadata = mcols(cds))
    message("Running PCA...")
    cds_wdoublets <- preprocess_cds(cds_wdoublets, num_dim = length(PCs), verbose = T, genes=ord_genes)
    pca.coord <- cds_wdoublets@reducedDims$PCA[ , PCs]
    cell.names <- rownames(colData(cds_wdoublets))
    nCells <- length(cell.names)
    rm(cds_wdoublets); gc() # Free up memory
    
    ## Compute PC distance matrix
    message("Calculating PC distance matrix...")
    dist.mat <- rdist(pca.coord)
    
    ## Compute pANN
    message("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      neighbor.names <- rownames(dist.mat)[neighbors]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
    }

    message("Classifying doublets..")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    colData(cds)[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(colData(cds)), 1]
    colData(cds)[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(cds)
}
  
  
  #' @importFrom pbmcapply pbmclapply
  #' @export
  #' 
  paramSweep_v3 <- function(cds, PCs=1:10, sct = FALSE, num.cores=detectCores()/2, genes=c("same", "all", "recalc"), ...) {
    ## Set pN-pK param sweep ranges
    pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.5,by=0.01))
    pN <- seq(0.05,0.3,by=0.05)
    
    ## Remove pK values with too few cells
    min.cells <- round(nrow(colData(cds))/(1-0.05) - nrow(colData(cds)))
    pK.test <- round(pK*min.cells)
    pK <- pK[which(pK.test >= 1)]
    
    ##get ordering genes
    og_done<-FALSE
    sweep.res.list = list()
    list.ind = 0
    dots <- list(...)
    sg_args <-c("logmean_ul", "logmean_ll", "top_n", "fit_min", "fit_max")
    cd_args<-c("id", "symbol_tag","method", "remove_outliers", "q")
    #gather relevant args
    rel_args<-dots[names(dots) %in% c(sg_args, cd_args)]
    if(!"top_n" %in% names(rel_args)){
      rel_args<-c(list(top_n=2000), rel_args)
    }
    if(length(genes)>1 & all(genes %in% c("same", "all", "recalc"))){
      genes <- "same"
    }
    if(length(genes)==1 & genes %in% c("same", "all", "recalc")){
      if(genes=="all") {
        message("Using all features")
        ord_genes <- rownames(fData(cds))
        og_done <-TRUE
      }
      if(genes=="recalc"){
        message("Recalculating ordering features using the following arguments:\nCalculate Dispersion:\n")
        print(rel_args[names(rel_args) %in% cd_args])
        message("\nSelect Genes:\n")
        print(rel_args[names(rel_args) %in% sg_args])
        rel_args<-c(list(cds=cds), rel_args)
        cds<-do.call(calculate_gene_dispersion, rel_args[names(rel_args) %in% c("cds", cd_args)])
        cds<-do.call(select_genes, rel_args[names(rel_args) %in% c("cds", sg_args)])
        ord_genes<-get_ordering_genes(cds)
        og_done <-TRUE
      }
      if(genes=="same"){
        message("Using existing ordering features")
        ord_genes<-get_ordering_genes(cds)
        og_done <-TRUE
      }
    }
    if(length(genes)>1 & all(genes %in% rownames(exprs(cds)))){
      ord_genes<-genes
      message("Using supplied ordering genes")
      og_done <-TRUE
    }
    if(!og_done & all(genes %in% rownames(exprs(cds)))){
      stop("Genes not found in cds; they must be rownames of exprs(cds)")
    }
    
    ## Down-sample cells to 10000 (when applicable) for computational effiency
    if (nrow(colData(cds)) > 10000) {
      real.cells <- rownames(colData(cds))[sample(1:nrow(colData(cds)), 10000, replace=FALSE)]
      data <- exprs(cds)[ , real.cells]
      n.real.cells <- ncol(data)
    }
    
    if (nrow(colData(cds)) <= 10000){
      real.cells <- rownames(colData(cds))
      data <- exprs(cds)
      n.real.cells <- ncol(data)
    }
    
    data<-data[rownames(data) %in% ord_genes,]
    ## Iterate through pN, computing pANN vectors at varying pK
    #no_cores <- detectCores()-1
    if(num.cores>1){
      #require(parallel)
      #cl <- makeCluster(num.cores)
      output2 <- pbmclapply(as.list(1:length(pN)),
                          FUN = parallel_paramSweep_v3,
                          n.real.cells,
                          real.cells,
                          pK,
                          pN,
                          data,
                          PCs,
                          mc.cores=num.cores)
      #stopCluster(cl)
    }else{
      output2 <- lapply(as.list(1:length(pN)),
                        FUN = parallel_paramSweep_v3,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        PCs)
    }
    
    ## Write parallelized output into list
    sweep.res.list <- list()
    list.ind <- 0
    for(i in 1:length(output2)){
      for(j in 1:length(output2[[i]])){
        list.ind <- list.ind + 1
        sweep.res.list[[list.ind]] <- output2[[i]][[j]]
      }
    }
    
    ## Assign names to list of results
    name.vec <- NULL
    for (j in 1:length(pN)) {
      name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
    }
    names(sweep.res.list) <- name.vec
    return(sweep.res.list)
    
  }

#' @importFrom fields rdist
#' @importFrom S4Vectors DataFrame
#' @export
#' 
parallel_paramSweep_v3 <- function(n, n.real.cells, real.cells, pK, pN, data, PCs)  {
  sweep.res.list = list()
  list.ind = 0
  ## Make merged real-artifical data
  message(paste("Creating artificial doublets for pN = ", pN[n]*100,"%",sep=""))
  n_doublets <- round(n.real.cells/(1 - pN[n]) - n.real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)
  
  ## Pre-process cds object
  message("Creating Monocle3 object with doublets...")
  cds_wdoublets <- new_cell_data_set(data_wdoublets, gene_metadata = DataFrame(gene_short_name=rownames(data), row.names = rownames(data)))
  message("Running PCA...")
  cds_wdoublets <- preprocess_cds(cds_wdoublets, num_dim = length(PCs), verbose = T)
  pca.coord <- cds_wdoublets@reducedDims$PCA[ , PCs]
  cell.names <- rownames(colData(cds_wdoublets))
  nCells <- length(cell.names)
  message("Calculating PC distance matrix...")
  nCells <- nrow(colData(cds_wdoublets))
  pca.coord <- cds_wdoublets@reducedDims$PCA[ , PCs]
  rm(cds_wdoublets)
  gc()
  dist.mat <- rdist(pca.coord)[,1:n.real.cells]
  
  ## Pre-order PC distance matrix prior to iterating across pK for pANN computations
  message("Defining neighborhoods...")
  for (i in 1:n.real.cells) {
    dist.mat[,i] <- order(dist.mat[,i])
  }
  
  ## Trim PC distance matrix for faster manipulations
  ind <- round(nCells * max(pK))+5
  dist.mat <- dist.mat[1:ind, ]
  
  ## Compute pANN across pK sweep
  message("Computing pANN across all pK...")
  for (k in 1:length(pK)) {
    print(paste("pK = ", pK[k], "...", sep = ""))
    pk.temp <- round(nCells * pK[k])
    pANN <- as.data.frame(matrix(0L, nrow = n.real.cells, ncol = 1))
    colnames(pANN) <- "pANN"
    rownames(pANN) <- real.cells
    list.ind <- list.ind + 1
    
    for (i in 1:n.real.cells) {
      neighbors <- dist.mat[2:(pk.temp + 1),i]
      pANN$pANN[i] <- length(which(neighbors > n.real.cells))/pk.temp
    }
    
    sweep.res.list[[list.ind]] <- pANN
    
  }
  return(sweep.res.list)
}


#' @importFrom ROCR prediction
#' @importFrom ROCR performance
#' @importFrom modes bimodality_coefficient
#' @importFrom KernSmooth bkde
#' @export
#' 
summarizeSweep <- function(sweep.list, GT = FALSE, GT.calls = NULL) {
  #require(KernSmooth); require(ROCR); require(modes)
  ## Set pN-pK param sweep ranges
  name.vec <- names(sweep.list)
  name.vec <- unlist(strsplit(name.vec, split="pN_"))
  name.vec <- name.vec[seq(2, length(name.vec), by=2)]
  name.vec <- unlist(strsplit(name.vec, split="_pK_"))
  pN <- as.numeric(unique(name.vec[seq(1, length(name.vec), by=2)]))
  pK <- as.numeric(unique(name.vec[seq(2, length(name.vec), by=2)]))
  
  ## Initialize data structure w/ or w/o AUC column, depending on whether ground-truth doublet classifications are available
  if (GT == TRUE) {
    sweep.stats <- as.data.frame(matrix(0L, nrow=length(sweep.list), ncol=4))
    colnames(sweep.stats) <- c("pN","pK","AUC","BCreal")
    sweep.stats$pN <- factor(rep(pN, each=length(pK), levels = pN))
    sweep.stats$pK <- factor(rep(pK, length(pN),levels = pK))
  }
  
  if (GT == FALSE) {
    sweep.stats <- as.data.frame(matrix(0L, nrow=length(sweep.list), ncol=3))
    colnames(sweep.stats) <- c("pN","pK","BCreal")
    sweep.stats$pN <- factor(rep(pN, each=length(pK), levels = pN))
    sweep.stats$pK <- factor(rep(pK, length(pN),levels = pK))
  }
  
  ## Perform pN-pK parameter sweep summary
  for (i in 1:length(sweep.list)) {
    res.temp <- sweep.list[[i]]
    
    ## Use gaussian kernel density estimation of pANN vector to compute bimodality coefficient
    gkde <- approxfun(bkde(res.temp$pANN, kernel="normal"))
    x <- seq(from=min(res.temp$pANN), to=max(res.temp$pANN), length.out=nrow(res.temp))
    sweep.stats$BCreal[i] <- bimodality_coefficient(gkde(x))
    
    if (GT == FALSE) { next }
    
    ## If ground-truth doublet classifications are available, perform ROC analysis on logistic
    ## regression model trained using pANN vector
    meta <- as.data.frame(matrix(0L, nrow=nrow(res.temp), ncol=2))
    meta[,1] <- GT.calls
    meta[,2] <- res.temp$pANN
    train.ind <- sample(1:nrow(meta), round(nrow(meta)/2), replace=FALSE)
    test.ind <- (1:nrow(meta))[-train.ind]
    colnames(meta) <- c("SinDub","pANN")
    meta$SinDub <- factor(meta$SinDub, levels = c("Doublet","Singlet"))
    model.lm <- glm(SinDub ~ pANN, family="binomial"(link='logit'), data=meta, subset=train.ind)
    prob <- predict(model.lm, newdata=meta[test.ind, ], type="response")
    ROCpred <- prediction(predictions=prob, labels=meta$SinDub[test.ind])
    perf.auc <- performance(ROCpred, measure="auc")
    sweep.stats$AUC[i] <- perf.auc@y.values[[1]]
  }
  
  return(sweep.stats)
  
}


#' @export
find.pK <- function(sweep.stats) {
  
  ## Implementation for data without ground-truth doublet classifications 
  '%ni%' <- Negate('%in%')
  if ("AUC" %ni% colnames(sweep.stats) == TRUE) {
    ## Initialize data structure for results storage
    bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=5))
    colnames(bc.mvn) <- c("ParamID","pK","MeanBC","VarBC","BCmetric")
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)
    
    ## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, "BCreal"])^2)
    }
    
    ## Plot for visual validation of BCmvn distribution
    # par(mar=rep(1,4))
    # x <- plot(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, pch=16, col="#41b6c4", cex=0.75)
    # x <- lines(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, col="#41b6c4")
    # print(x)
    # 
    return(bc.mvn)
    
  }
  
  ## Implementation for data with ground-truth doublet classifications (e.g., MULTI-seq, CellHashing, Demuxlet, etc.)
  if ("AUC" %in% colnames(sweep.stats) == TRUE) {
    ## Initialize data structure for results storage
    bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=6))
    colnames(bc.mvn) <- c("ParamID","pK","MeanAUC","MeanBC","VarBC","BCmetric")
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)
    
    ## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanAUC[x] <- mean(sweep.stats[ind, "AUC"])
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, "BCreal"])^2)
    }
    
    ## Plot for visual validation of BCmvn distribution
    # par(mar=rep(1,4))
    # x <- plot(x=bc.mvn$ParamID, y=bc.mvn$MeanAUC, pch=18, col="black", cex=0.75,xlab=NA, ylab = NA)
    # x <- lines(x=bc.mvn$ParamID, y=bc.mvn$MeanAUC, col="black", lty=2)
    # par(new=TRUE)
    # x <- plot(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, pch=16, col="#41b6c4", cex=0.75)
    # axis(side=4)
    # x <- lines(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, col="#41b6c4")
    # print(x)
    
    return(bc.mvn)
    
  }
}

#' @export
modelHomotypic <- function(annotations) {
  anno.freq <- table(annotations)/length(annotations)
  x <- sum(anno.freq^2)
  return(x)
}

seurat_paramSweep_v3<-function (seu, PCs = 1:10, sct = FALSE, num.cores = 1) 
{
  require(Seurat)
  require(fields)
  pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.5, by = 0.05)
  min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]
  orig.commands <- seu@commands
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 
                                                 10000, replace = FALSE)]
    data <- seu@assays$RNA@counts[, real.cells]
    n.real.cells <- ncol(data)
  }
  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }
  if (num.cores > 1) {
    require(parallel)
    #cl <- makeCluster(num.cores)
    output2 <- pbmcapply::pbmclapply(as.list(1:length(pN)), FUN = DoubletFinder::parallel_paramSweep_v3, 
                        n.real.cells, real.cells, pK, pN, data, orig.commands, 
                        PCs, sct, mc.cores = num.cores)
    #stopCluster(cl)
  }
  else {
    output2 <- lapply(as.list(1:length(pN)), FUN = DoubletFinder::parallel_paramSweep_v3, 
                      n.real.cells, real.cells, pK, pN, data, orig.commands, 
                      PCs, sct)
  }
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, 
                                  sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}
