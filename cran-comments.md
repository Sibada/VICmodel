## Resubmission
This is a resubmission. In this version I have:

* Fix the significant warnings of the compilation of VIC source codes under
win-builder environment.


## Test environments
* Mac OS X 10.13.3 (on travis-ci), R 3.5.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.0
* win-builder (devel and release)

## R CMD check results

There was 2 NOTE for ubuntu environment:

* This is a new release.

* checking installed package size ... NOTE
  installed size is  5.3Mb
  sub-directories of 1Mb or more:
    libs   4.7Mb

  This is because the VIC model is a complex model and the core has large code. 
  That is unavoidable but have little impact for users, therefore we suggest that
  it is not a significant problem.

There was 1 NOTE for Mac OS X environment:

* This is a new release.

There was 3 NOTE for win-builder environment:

* This is a new release.

* running examples for arch 'i386' ... [14s] NOTE
Examples with CPU or elapsed time > 10s
    user system elapsed
vic 12.9   0.07   13.01

* running examples for arch 'x64' ... [12s] NOTE
Examples with CPU or elapsed time > 10s
     user system elapsed
vic 11.03   0.03   11.09
  
  This might be mainly because of the differences between the C source codes of VIC running on the linux and windows   platform. Moreover, the VIC model is a complex landsurface scheme and hence would require more computation than the   other landsurface models. We suggest that it is normal thus not a significant problem.
  
## Reverse dependencies

This is a new release, so there are no reverse dependencies.
