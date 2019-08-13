library(monocle3)
library(m3addon)
reticulate::py_config()
library(magrittr)
library(reticulate)
library(Rcpp)
roxygen2::roxygenize(".")
usethis::use_build_ignore("debug.R")

cds<-readRDS("data/m3cds.RDS")


sm<-as.matrix(normalized_counts(cds[1:1000,sample(ncol(cds), 1000)]))
dim(sm)
x <- matrix(rnorm(100), nrow = 20)
coords<-spRing(x, method="euclidean")

bs = backspin(sm)

library(RcppCNPy)
npySave("data/sm.npy", sm)
getwd()

#in python
import numpy as np
sm = np.load('/Users/sfurla/Box Sync/PI_FurlanS/computation/Rproj/m3addon/data/sm.npy')
backSPIN(sm)
#seems to work


# system.time({
# rs<-m3addon:::rowStdDev(exprs(cds))
# })
# 
# system.time({
# cs<-m3addon:::colStdDev(t(exprs(cds)))
# })
# all(rs[,1]==cs[1,])

system.time({
  rt<-m3addon:::rowStdDev(exprs(cds))
})

system.time({
  rr<-rowSds(as.matrix(exprs(cds)))
})
# all(as.numeric(rt[1,])[1:20]==rr[1:20])


Rcpp::sourceCpp("src/scores.cpp")

cds<-calculate_gene_dispersion(cds, q=5)
plot_gene_dispersion(cds)
cds<-select_genes(cds)
plot_gene_dispersion(cds)

ord_genes<-get_ordering_genes(cds)
cds<-preprocess_cds(cds, use_genes = ord_genes, verbose = T, num_dim = 100)
plot_pc_variance_explained(cds)
cds<-reduce_dimension(cds, reduction_method = "UMAP", num_dim = 35, verbose=T, cores = detectCores())
plot_cells(cds, color_cells_by = "Group")


levels(factor(pData(cds)$Group))
cdsDE<-cds[,pData(cds)$Group %in% c("Colon_Donor", "Colon_Host")]
levels(factor(pData(cdsDE)$Group))

gene_fits <-fit_models(cdsDE[1:10,], model_formula_str = "~Group", verbose = TRUE, cores = detectCores())


###########AF#############
s <- matrix(seq(0, 100, by = .0001), ncol = 1)
rbenchmark::benchmark(Arma = put_option_pricer_arma(s, 60, .01, .02, 1, .05),
                      AF = put_option_pricer_af(s, 60, .01, .02, 1, .05),
                      order = "relative", 
                      replications = 100)[,1:4]

Rcpp::sourceCpp("src/armatut.cpp")
Rcpp::sourceCpp("src/aftut.cpp", rebuild = T)
