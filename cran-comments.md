# VICmodel 0.1.1

* Fix some codes with bug or potential problem of the VIC source codes, including:
  
Comment out functions like backtrace/backtrace_symbols.

Macro definition `ERROR` is renamed as `VIC_ERROR`.

Comment out the macros that are not used in this package such as those of openmp.

Comment out the unnecessary print functions.

* Fix the data error of the `STEHE` sample data for routing model.

* Add the observed streamflow data to `STEHE` sample data.

## Test environments
* Mac OS X 10.13.3 (on travis-ci), R 3.5.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.0
* win-builder (devel and release)

## R CMD check results

There was 1 NOTE for ubuntu environment:

* This is a new release.

There was 1 NOTE for Mac OS X environment:

* This is a new release.

There was 1 NOTE for win-builder environment:

* This is a new release.

  
## Reverse dependencies

This is a new release, so there are no reverse dependencies.
