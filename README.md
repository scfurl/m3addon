# m3addon

This package adds to the popular "monocle3" for use with the R computing environment

## Getting Started

These instructions will get you a copy of the project up and running on your local machine.

### Prerequisites

Imports: 
  Rcpp (>= 1.0.0),
  monocle3 (>= 0.1.0),
  ggplot2 (>= 3.0.0),
  h5 (>= 0.9.9)

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

Check out the following to get Rcpp up and running on your system:

https://github.com/RcppCore/Rcpp/issues/922

## Using this software

### Load h5

A monocle3 cds can be created from a vector of folders as such:

```
cds<-load_cellranger_data_h5(folders)

```

Data generated using the cellranger 'aggregate' function can be used to create a monocle3 cds as follows:  (note that only 1 folder is currently supported)

```
cds<-load_cellranger_data_h5(folder, aggregated=T)

```

## Authors

* **Scott Furlan** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used - especially Trapnell lab folks
