## Resubmission
This is a resubmission. In this version I have:

* Added small executable examples in my Rd-files to show the use of the exported function.
* Added a citation in the DESCRIPTION file to an arXiv paper describing the methods in my package.

## Test environments
* ubuntu 16.04.6 LTS (on travis-ci) with R 4.02
* win-builder (checked with both devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* check for GNU extensions in Makefiles ... 
  GNU make is a SystemRequirements.

  My package depends on RcppParallel which depends on GNU make extensions.
  This is a note that is only thrown when running the check on Linux, and is
  not shown with the win-builder check. To my knowledge, this is a standard
  note to expect for a package that depends on RcppParallel.

## Downstream dependencies
There are currently no downstream dependencies for this package.
