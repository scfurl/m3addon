# m3addon

This package adds to the popular "monocle3" package for use in the R computing environment.  The goal of this repo is to coalesce working R implementations of popular functions used in single cell RNA seq analysis but are inaccessible to the R user because they are split across programming languages or are not particularly well suited for working in an interactive environment. This package includes R implementations for functions from other packages such as spRing plots, backSPIN, scrublet.  The citations associated with these functions are included in the documentation.  It also adds some functionality to monocle3 including tools that can read 10X data from h5 files and import them directly into a monocle3 cell_data_set.  See below.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine.

### Prerequisites

FYI, imports and requires: 
  Rcpp (>= 1.0.0),
  monocle3 (>= 0.1.0),
  ggplot2 (>= 3.0.0),
  h5 (>= 0.9.9),
  reticulate (>=1.12)

### Installing

Check out the following to install monocle3:  https://cole-trapnell-lab.github.io/monocle3/monocle3_docs/#installing-monocle

Once you have it up and running, perform the following

```
devtools::install_github('scfurl/m3addon')
```

If you are using macOS get the following error: 

```
/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/math.h:301:15: fatal error: 'math.h' file not found
```

Check out the following to get Rcpp compilation up and running on your macOS:

https://github.com/RcppCore/Rcpp/issues/922

## Using this software

### Load h5

A monocle3 cds can be created from a vector of cellranger folders as such:

```
cds<-load_cellranger_data_h5(folders)

```

Data generated using the cellranger 'aggregate' function can be used to create a monocle3 cds as follows:  (note that only 1 folder is currently supported)

```
cds<-load_cellranger_data_h5(folder, aggregated=T)

```

The resulting cds will have a colData object that incorporates information from the cellranger output file 'aggregate.csv'.  As such, this csv file as named should be included in the folder.  Further details for this function and all others can be found in the package documentation.

## Authors

* **Scott Furlan** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used - especially Trapnell lab folks
