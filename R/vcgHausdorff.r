#' calculate the Hausdorff beteween 2 meshes
#'
#' calculate the Hausdorff beteween two triangular meshes.
#'
#' @param ref triangular mesh (object of class 'mesh3d').
#' @param target triangular mesh (object of class 'mesh3d').
#' @param nSamples set the required number of samples.
#' @param nSamplesArea set the required number of samples per area unit, override nSamples.
#' @param vertSamp logical: if TRUE, disable vertex sampling.
#' @param edgeSamp logical: if TRUE, disable edge sampling.
#' @param faceSamp logical: if TRUE, disable face sampling.
#' @param unrefVert logical: if FLASE, ignore unreferred vertices.
#' @param samplingType SS (similar triangles sampling), SD (subdivision sampling), MC (montecarlo sampling).
#' @param searchStruct SGIRD (static Uniform Grid), OCTREE, AABB (AxisAligned Bounding Box Tree), HGRID (Hashed Uniform Grid).
#' @param from numeric: minimum value for color mapping.
#' @param to numeric: maximum value for color mapping.
#'
#' @return
#' \item{ForwardSampling, BackwardSampling}{lists containing information about forward sampling with the following entries}
#' \itemize{
#' \item{maxdist}{maximal Hausdorff distance}
#' \item{meandist}{mean Hausdorff distance}
#' \item{RMSdist}{RMS of the Hausdorff distances}
#' \item{nvbsamples}{number of vertices sampled}
#' \item{nsamples}{number of samples}
#' }
#' \item{mesh1, mesh2}{meshes with color coded distances and an additional entry called quality containing the sampled per-vertex distances}
#' \item{forward_hist, backward_hist}{Matrices tracking the sampling results}
#'
#' @export vcgHausdorff
# vcgHausdorff <- function(ref, target) {
vcgHausdorff <- function(ref, target, nSamples=0, nSamplesArea=0, vertSamp=F, edgeSamp=F, faceSamp=F, unrefVert=T,samplingType=c("SS","MC","SD"),searchStruct=c("SGRID","AABB","OCTREE","HGRID"), from=0, to=0){
	if (!inherits(ref,"mesh3d"))
		stop("argument 'ref' needs to be object of class 'mesh3d'")
	if (!inherits(target,"mesh3d"))
		stop("argument 'target' needs to be object of class 'mesh3d'")

	samplingType <- match.arg(samplingType[1],c("SS","MC","SD") )
	if (samplingType == "MC")
		samplingType <- 0
	if (samplingType == "SD")
		samplingType <- 1
	if (samplingType == "SS")
		samplingType <- 2

	searchStruct <- match.arg(searchStruct[1],c( "SGRID", "AABB", "OCTREE", "HGRID") )
	if (searchStruct == "AABB")
		searchStruct <- 0
	if (searchStruct == "SGRID")
		searchStruct <- 1
	if (searchStruct == "HGRID")
		searchStruct <- 2
	if (searchStruct == "OCTREE")
		searchStruct <- 3

   tmp <- .Call("Rhausdorff",ref$vb[1:3,],ref$it-1L,target$vb[1:3,],target$it-1L, vertSamp, edgeSamp, faceSamp, unrefVert, samplingType, nSamples, nSamplesArea, from, to, searchStruct)
        tmp$mesh1 <- updateColorFromVector(tmp$mesh1$mesh,tmp$mesh1$colors)
        tmp$mesh2 <- updateColorFromVector(tmp$mesh2$mesh,tmp$mesh2$colors)
        
        
   return(tmp)
}


updateColorFromVector <- function(out,colors) {
    out$material <- list()
            if (length(colors)) {
                colvec <- matrix(colors,3,(length(colors)/3))
                colvec <- rgb(colvec[1,],colvec[2,],colvec[3,],maxColorValue=255)
                colfun <- function(x)
                    {
                        x <- colvec[x]
                        return(x)
                    }
                out$material$color <- matrix(colfun(out$it),dim(out$it))
            }
    class(out) <- "mesh3d"
            
    return(out)
}
