% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/publication_plotting.R
\name{red_dim_plot}
\alias{red_dim_plot}
\title{Reduced Dimensionality Plot}
\usage{
red_dim_plot(
  plot,
  method = "UMAP",
  size = 20,
  ICON_DIR = system.file("icons", package = "m3addon"),
  correction_x = 0,
  correction_y = 0
)
}
\arguments{
\item{plot}{A ggplot2 object}

\item{method}{method indicating which icon to use.  Currently 'tSNE', 'PCA', and 'UMAP' supported}

\item{correction_x, }{integer to offset icon in x axis}

\item{correction_y, }{integer to offset icon in y axis}

\item{widths}{A vector of length 2 specifying the size of the plot grid and legend grid, respectively.}
}
\value{
grid arranged plot using split_plot function
}
\description{
Adds an icon to the lower left hand corner of a plot indicating reducted dimensionality method used
}
\references{
Puram, S. V. et al. Single-Cell Transcriptomic Analysis of Primary and 
Metastatic Tumor Ecosystems in Head and Neck Cancer. Cell 171, 1611.e1–1611.e24 (2017).
}
