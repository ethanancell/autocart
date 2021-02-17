## Resubmission
This is a resubmission. In this version I have:

* Added small executable examples in my Rd-files to show the use of the exported function.
* Added a citation in the DESCRIPTION file to an arXiv paper describing the methods in my package.

## Test environments
* local Linux install (Pop!\_OS 20.10), R 4.0.2
* win-builder (checked with both devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 4 NOTEs:

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
  
* running examples for arch 'i386' ... [384s] NOTE
  Examples with CPU (user + system) or elapsed time > 10s
  
  The bulk of this long example run time comes from the "autotune" (tuning hyperparameters)
  and "autoforest" (ensemble of autocart trees) functions. In each of these functions,
  many autocart trees are created which takes a lot of computation. To get my runtime
  below 10 seconds would require me to remove the code example entirely or strip down
  the code example so much that it would be difficult to learn how to apply my hyperparameter
  tuning function in practice from the code example. This is also a NOTE that is not
  thrown when I check my package on Linux, as my run times seem to be fine when running locally.
  
* running examples for arch 'x64' ... [407s] NOTE
  Examples with CPU (user + system) or elapsed time > 10s
  
  My explanation for this note is the same as the note above with arch i386. Again,
  this is not a NOTE that is thrown when I check the package on Linux, as my local
  machine seems to run the example code quick enough.

## Downstream dependencies
There are currently no downstream dependencies for this package.
