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
#' @param saveMesh logical: if TRUE, save a mesh with error as per-vertex colour and quality.
#' @param from: numeric: minimum value for color mapping.
#' @param to: numeric: maximum value for color mapping.
#' @param writeHist logical: if TRUE, write files with histograms of error distribution.
#'
#' @return
#' \item{maxdist}{maximal Hausdorff distance}
#' \item{meandist}{mean Hausdorff distance}
#' \item{RMSdist}{RMS of the Hausdorff distances}
#' \item{nvbsamples}{number of vertices sampled}
#' \item{nsamples}{number of samples}
#'
#' @export vcgHausdorff
# vcgHausdorff <- function(ref, target) {
vcgHausdorff <- function(ref, target, nSamples=0, nSamplesArea=0, vertSamp=F, edgeSamp=F, faceSamp=F, unrefVert=T,samplingType=c("SS","MC","SD"),searchStruct=c("SGRID","AABB","OCTREE","HGRID"), saveMesh=F, from=0, to=0, writeHist=F){
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

   tmp <- .Call("Rhausdorff",ref$vb[1:3,],ref$it-1L,target$vb[1:3,],target$it-1L, vertSamp, edgeSamp, faceSamp, unrefVert, samplingType, nSamples, nSamplesArea, saveMesh, from, to, writeHist, searchStruct)
   return(tmp)
}