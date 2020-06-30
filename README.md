# Autocart
A modified regression tree R package so that spatial autocorrelation can be tracked at every step of the splitting.

Existing classification and regression tree (CART) R packages allow for custom splitting functions, but not at high enough of a level to
where geographical coordinates can be tracked at every stage to include custom splitting that accounts for spatial autocorrelation.

This package is an implementation of the CART algorithm that allows for tracking geographical location in the observations at every
stage of the tree creation.

## Installation from source
To install this package from the source code, make sure that you have the devtools package downloaded by using `install.packages("devtools")`.
### Windows
You must have Rtools downloaded so that the C++ source can be compiled. The most recent version of rtools can be found [here](https://cran.r-project.org/bin/windows/Rtools/)
### macOS
Install the Xcode command line tools with `xcode-select --install` in the shell. You may need to register as an Apple developer first.
### Linux
To compile the C++ code, you must also have the R development tools, which can be installed by installing the r-base-dev package.

### After downloading compiler
To install this package in an R environment, use `devtools::install_github("ethanancell/autocart")`

## Usage
To get started after installation, view the introductory autocart vignette by using `vignette("autocart-intro")`

## License
[MIT](LICENSE)
