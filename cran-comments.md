## Resubmission
This is a resubmission. In this version I have:

* Replaced the "\dontrun{}" by "\donttest{}" in the Rd files.

* Refered the copyright holders in the inst/COPYRIGHTS, including:
  
  Copyright (C) 2016 The Computational Hydrology Group, Department of 
  Civil and Environmental Engineering, University of Washington.
  
  Copyright (C) 2014 The Land Surface Hydrology Group, Department of
  Civil and Environmental Engineering, University of Washington.
  (*Note*: this is the same organization of the Computational Hydrology
  Group, for the former renamed as the latter, therefore we only retain
  those of Computational Hydrology Group of UW)
  
  Copyright (c) 2010, Zed A. Shaw and Mongrel2 Project Contributors.
  
* Add some relevant contributers as author marked as "ctb".
  
* Add some necessary references.

Thanks for your review.

Best regards

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

There was 1 NOTE for win-builder environment:

* This is a new release.

  
## Reverse dependencies

This is a new release, so there are no reverse dependencies.
