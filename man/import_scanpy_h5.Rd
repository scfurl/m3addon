% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scanpy.R
\name{import_scanpy_h5}
\alias{import_scanpy_h5}
\title{Imports a scanpy h5 object}
\usage{
import_scanpy_h5(
  file,
  raw = TRUE,
  rowname_col = "GeneID",
  gene_short_name_col = "GeneName",
  colors = NULL,
  PCA = T,
  UMAP = TRUE,
  neighbors = T,
  var_features = T
)
}
\arguments{
\item{file}{h5 scanpy object}

\item{raw}{data contained in the raw slot - See scanpy documentation for more info}

\item{rowname_col}{feature column for rownames (should contain unique values - suggest ENSID)}

\item{gene_short_name_col}{feature column for gene_short_name (suggest Symbol)}

\item{colors}{name or names of color palatte to import}

\item{exprs_bin}{Boolean whether to bin genes by mean expression}

\item{exprs_cuts}{Numeic indicating number of bins if using exprs_bin}
}
\value{
Updated cell_data_set object
}
\description{
read in a Scanpy h5 object.
}
