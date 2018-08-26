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

There was 2 NOTE for Mac OS X environment:

* This is a new release.

* checking whether package ¡®VICmodel¡¯ can be installed ... NOTE
  Found the following warnings:
    vic/vic_run/src/CalcBlowingSnow.c:308:1: warning: control may reach end of non-void function [-Wreturn-type]
    vic/vic_run/src/CalcBlowingSnow.c:487:1: warning: control may reach end of non-void function [-Wreturn-type]
    vic/drivers/shared_all/src/input_tools.c:44:1: warning: control may reach end of non-void function [-Wreturn-type]
    vic/drivers/shared_all/src/input_tools.c:135:1: warning: control may reach end of non-void function [-Wreturn-type]
    vic/drivers/shared_all/src/input_tools.c:169:1: warning: control may reach end of non-void function [-Wreturn-type]
    vic/drivers/shared_all/src/input_tools.c:224:1: warning: control may reach end of non-void function [-Wreturn-type]
    vic/drivers/shared_all/src/input_tools.c:276:1: warning: control may reach end of non-void function [-Wreturn-type]
    vic/drivers/shared_all/src/input_tools.c:300:1: warning: control may reach end of non-void function [-Wreturn-type]

  They are originated from the [source code of the VIC model](https://github.com/UW-Hydro/VIC). Although has those
  warnings when compiling as shown above, those codes has been widely applied and verified for several applications for a long time e.g. scienctific research, development of disaster monitoring system and land data assimilation system (LDAS) since the VIC model is firstly developed in 1994 by [Liang et al. (1994)](dx.doi.org/10.1029/94JD00483). Therefore, we suggest this is not a significant problem.

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
  
  This might be mainly because of the differences between the C source codes of VIC running on the linux and windows platform. Moreover, the VIC model is a complex landsurface scheme and hence would require more computation than the other landsurface models. We suggest that it is normal thus not a significant problem.
  
## Reverse dependencies

This is a new release, so there are no reverse dependencies.
