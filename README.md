
Rvcg [![Unit Tests](https://github.com/zarquon42b/Rvcg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zarquon42b/Rvcg/actions/workflows/R-CMD-check.yaml)  [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)  [![Downloads](https://cranlogs.r-pkg.org/badges/last-month/Rvcg?color=brightgreen)](https://cranlogs.r-pkg.org/badges/last-month/Rvcg)
====
__Rvcg__ is an R-package providing methods for manipulations on triangular meshes by using the API of the [VCGLIB](http://vcg.sf.net/) library.

#### Installation of the R-package Rvcg from CRAN: ####

Within R:
       
       install.packages("Rvcg")


#### Installation of the-R package "Rvcg" (latest development code) using *devtools*:: ####

##### install prerequisites #####

1. install *devtools* from within R (Ubuntu/Debian users will have to install *libcurl4-gnutls-dev* beforehand):

        install.packages("devtools")

2. Install build environment
    * **Windows:** Install latest version of *[Rtools](http://cran.r-project.org/bin/windows/Rtools)*
During installation of *Rtools* make sure to install the *toolchain*, and to select *"Edit the system path"* (and confirming the installers suggestions).
    * **OSX:** Install *[XCODE](https://developer.apple.com/xcode/)*

##### install Rvcg #####
Run the following command in R:
        
        require(devtools)
        install_github("zarquon42b/Rvcg", local=FALSE)
    
### Conda Version 
Rvcg has been ported to the python packaging system Conda thanks to [@dipterix](https://github.com/username@dipterix) and can be found here [https://anaconda.org/conda-forge/r-rvcg](https://anaconda.org/conda-forge/r-rvcg)

From an active conda enironment, it can be installed with the command
 
 ```conda install conda-forge::r-rvcg ```
