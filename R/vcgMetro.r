#' evaluate the difference between two triangular meshes.
#'
#' Implementation of the command line tool "metro" to evaluate the difference between two triangular meshes.
#'
#' @param mesh1 triangular mesh (object of class 'mesh3d').
#' @param mesh2 triangular mesh (object of class 'mesh3d').
#' @param nSamples set the required number of samples if 0, this will be set to approx. 10x the face number.
#' @param nSamplesArea set the required number of samples per area unit, override nSamples.
#' @param vertSamp logical: if FALSE, disable vertex sampling.
#' @param edgeSamp logical: if FALSE, disable edge sampling.
#' @param faceSamp logical: if FALSE, disable face sampling.
#' @param unrefVert logical: if FALSE, ignore unreferred vertices.
#' @param samplingType set the face sampling mode. options are: SS (similar triangles sampling), SD (subdivision sampling), MC (montecarlo sampling).
#' @param searchStruct set search structures to use. options are: SGIRD (static Uniform Grid), OCTREE, AABB (AxisAligned Bounding Box Tree), HGRID (Hashed Uniform Grid).
#' @param from numeric: minimum value for color mapping.
#' @param to numeric: maximum value for color mapping.
#' @param colormeshes if TRUE, meshes with vertices colored according to distance are returned
#' @param silent logical: if TRUE, output to console is suppressed.
#'
#' @return
#' \item{ForwardSampling, BackwardSampling}{lists containing information about forward (mesh1 to mesh2) and backward (mesh2 to mesh1) sampling with the following entries}
#' \itemize{
#' \item{ \code{maxdist} maximal Hausdorff distance}
#' \item{ \code{meandist} mean Hausdorff distance}
#' \item{ \code{RMSdist} RMS of the Hausdorff distances}
#' \item{ \code{area} mesh area (of \code{mesh1} in \code{ForwardSampling} and  \code{mesh2} in  \code{BackwardSampling})}
#' \item{ \code{RMSdist} RMS of the Hausdorff distances}
#' \item{ \code{nvbsamples} number of vertices sampled}
#' \item{ \code{nsamples} number of samples}
#' }
#' \item{distances1, distances2}{vectors containing vertex distances from mesh1 to mesh2 and mesh2 to mesh1.}
#' \item{forward_hist, backward_hist}{Matrices tracking the sampling results}
#' if colormeshes == TRUE
#' \item{mesh1, mesh2}{meshes with color coded distances and an additional entry called quality containing the sampled per-vertex distances}
#' @examples
#' require(Morpho)
#' data(humface)
#' data(dummyhead)
#' ## align humface to dummyhead.mesh
#' humfalign <- rotmesh.onto(humface,humface.lm,dummyhead.lm)
#' samp <- vcgMetro(humfalign$mesh,dummyhead.mesh,faceSamp=FALSE,edgeSamp=FALSE)
#' ## create heatmap using Morpho's meshDist function
#' \dontrun{
#' ## create custom heatmaps based on distances
#' mD <- meshDist(humfalign$mesh,distvec=samp$distances1)
#' }
#' @note this is a straightforward implementation of the command line tool metro \url{http://vcg.isti.cnr.it/vcglib/metro.html}
#' @references
#' P. Cignoni, C. Rocchini and R. Scopigno. Metro: measuring error on simplified surfaces. Computer Graphics Forum, Blackwell Publishers, vol. 17(2), June 1998, pp 167-174
#' @export
vcgMetro <- function(mesh1, mesh2, nSamples=0, nSamplesArea=0, vertSamp=TRUE, edgeSamp=TRUE, faceSamp=TRUE, unrefVert=FALSE,samplingType=c("SS","MC","SD"),searchStruct=c("SGRID","AABB","OCTREE","HGRID"), from=0, to=0, colormeshes=FALSE, silent=FALSE){
	if (!inherits(mesh1,"mesh3d"))
		stop("argument 'mesh1' needs to be object of class 'mesh3d'")
	if (!inherits(mesh2,"mesh3d"))
		stop("argument 'mesh2' needs to be object of class 'mesh3d'")

	samplingType <- match.arg(samplingType[1],c("SS","MC","SD","NS") )
	if (samplingType == "MC")
		samplingType <- 0
	if (samplingType == "SD")
		samplingType <- 1
	if (samplingType == "SS")
		samplingType <- 2
        ## if (samplingType == "NS")
	## 	samplingType <- 3

	searchStruct <- match.arg(searchStruct[1],c( "SGRID", "AABB", "OCTREE", "HGRID") )
	if (searchStruct == "AABB")
		searchStruct <- 0
	if (searchStruct == "SGRID")
		searchStruct <- 1
	if (searchStruct == "HGRID")
		searchStruct <- 2
	if (searchStruct == "OCTREE")
		searchStruct <- 3

        tmp <- .Call("Rmetro",mesh1,mesh2, vertSamp, edgeSamp, faceSamp, unrefVert, samplingType, nSamples, nSamplesArea, from, to, searchStruct,colormeshes,silent)
        tmp$mesh1 <- updateColorFromVector(tmp$mesh1$mesh,tmp$mesh1$colors)
        tmp$mesh2 <- updateColorFromVector(tmp$mesh2$mesh,tmp$mesh2$colors)
        
        
   return(tmp)
}


updateColorFromVector <- function(out,colors) {
    out$material <- list()
            if (length(colors)) {
                colvec <- matrix(colors,3,(length(colors)/3))
                colvec <- rgb(colvec[1,],colvec[2,],colvec[3,],maxColorValue=255)
                #colfun <- function(x)
                #    {
                #        x <- colvec[x]
                #        return(x)
                #    }
                out$material$color <-colvec
                    #matrix(colfun(out$it),dim(out$it))
            }
    class(out) <- "mesh3d"
            
    return(out)
}
