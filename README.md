Rvcg
====
__Rvcg__ is an R-package providing methods for manipulations on triangular meshes by using the API of the [VCGLIB](http://vcg.sf.net/) library.

#### Installation of the R-package Rvcg from CRAN: ####

Within R:
       
       install.packages("Rvcg")

#### Installation of the R-package "Rvcg": ####
0. Make sure to work with the latest version of R and install dependencies (type the following commands into your R terminal):                
        
	    install.packages("RcppEigen")

1. Download the version suitable for your OS from [here](https://github.com/zarquon42b/Rvcg/releases/). Either the compiled package (for Windows and OS X) or the source tarball (Linux).

2. Installation command from within R: 
   
        install.packages("Path_to_downloaded_package_Rvcg[Version_OS]",repos=NULL)

3. check if the package can be loaded:
        
        load package: library(Rvcg)

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
    
