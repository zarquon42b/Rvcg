#' Interface between R and vcglib libraries for mesh operations
#' 
#' Provides meshing functionality from vcglib (meshlab) for R. E.g. mesh
#' smoothing, mesh decimation, closest point search.
#' 
#' \tabular{ll}{
#' Package: \tab Rvcg\cr
#' Type: \tab Package\cr
#' Version: \tab 0.4\cr
#' Date: \tab 2013-08-02\cr
#' License: \tab GPL\cr
#' LazyLoad: \tab yes\cr }
#' 
#' @name Rvcg-package
#' @aliases Rvcg-package Rvcg
#' @docType package
#' @author Stefan Schlager
#' 
#' Maintainer: Stefan Schlager <zarquon42@@gmail.com>
#' @references To be announced
#' @keywords package
#' @useDynLib Rvcg
NULL

#' Example mesh
#'
#' A triangular mesh representing a human face - called by \code{data(humface)}
#' 
#' @name humface
#' @docType data
#' 
NULL
#' landmarks on mesh 'humface'
#' @usage data(humface)
#' landmarks on mesh 'humface'- called by \code{data(humface)}
#' @name humface.lm
#' @docType data
#' 
NULL

#' dummyhead mesh - called by data(dummyhead)
#'
#' A triangular mesh representing a dummyhead - called by \code{data(dummyhead)}
#' 
#' @name dummyhead.mesh
#' @docType data
#' 
NULL
#' landmarks on mesh 'dummyhead'
#' 
#' landmarks on mesh 'dummyhead'- called by \code{data(dummyhead)}
#'
#' @name dummyhead.lm
#' @docType data
#' 
NULL
