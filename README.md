# m3addon

This package adds to the popular "monocle3" for use with the R computing environment

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

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


## Authors

* **Scott Furlan** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used - especially Trapnell lab folks
