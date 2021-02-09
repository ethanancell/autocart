## Resubmission
This is a resubmission. In this version I have:

* Added small executable examples in my Rd-files to show the use of the exported function.

## Test environments
* local Linux install (Pop!\_OS 20.10), R 4.0.2
* win-builder (checked with both devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs:

* check for GNU extensions in Makefiles ... 
  GNU make is a SystemRequirements.

  The package depends on RcppParallel which depends on GNU make extensions.
  This is a note that is only thrown when running the check on Linux.

* checking installed package size ...
  installed size is 8.0Mb
  sub-directories of 1Mb or more:
    libs  7.8Mb
    
  In the installed package file under "libs", there is only one file which is
  the shared library object "autocart.so", which is the file of concern. This
  package relies heavily on a lot of compiled code, so it makes sense that
  this file might be a bit larger. I do not know of a way to reduce this
  filesize other than removing some of the functionality of my package. This
  is also a note that is not thrown on the Windows check, this is only seen
  when running the check with Linux.

## Downstream dependencies
There are currently no downstream dependencies for this package.
