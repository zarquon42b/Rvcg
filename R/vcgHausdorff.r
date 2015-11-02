#' calculate the Hausdorff beteween 2 meshes
#'
#' calculate the Hausdorff beteween two triangular meshes.
#'
#' @param ref triangular mesh (object of class 'mesh3d')
#' @param target triangular mesh (object of class 'mesh3d')
#'
#' @return
#' \item{maxdist}{maximal Hausdorff distance}
#' \item{meandist}{mean Hausdorff distance}
#' \item{RMSdist}{RMS of the Hausdorff distances}
#' \item{nvbsamples}{number of vertices sampled}
#' \item{nsamples}{number of samples}
#'
#' @examples
#'
#' data(humface)
#' curv <- vcgCurve(humface)
#' ##visualise per vertex mean curvature
#' \dontrun{
#' require(Morpho)
#' meshDist(humface,distvec=curv$meanvb,from=-0.2,to=0.2,tol=0.01)
#' }
#' @export vcgHausdorff
vcgHausdorff <- function(ref, target)
	{
		if (!inherits(ref,"mesh3d"))
			stop("argument 'ref' needs to be object of class 'mesh3d'")
		if (!inherits(target,"mesh3d"))
			stop("argument 'target' needs to be object of class 'mesh3d'")
		meshes <- list(ref,target)
		vb <- list()
		it <- list()
		for (i in 1:length(meshes)) {
			meshes[[i]] <- meshintegrity(meshes[[i]],facecheck=TRUE)
			vb[[i]] <- meshes[[i]]$vb[1:3,,drop=FALSE]
			it[[i]] <- meshes[[i]]$it - 1
			storage.mode(it[[i]]) <- "integer"
		}

		# ref <- meshintegrity(ref,facecheck=TRUE)
		# target <- meshintegrity(target,facecheck=TRUE)
		# vb0 <- ref$vb0[1:3,,drop=FALSE]
		# it0 <- ref$it0 - 1
		# dimit0 <- dim(it0)[2]
		# dimvb0 <- dim(vb0)[2]
		# storage.mode(it0) <- "integer"

		# vb1 <- target$vb1[1:3,,drop=FALSE]
		# it1 <- target$it1 - 1
		# dimit1 <- dim(it1)[2]
		# dimvb1 <- dim(vb1)[2]
		# storage.mode(it) <- "integer"
		tmp <- .Call("Rhausdorff",vb[1],it[1],vb[2],it[2])
		return(tmp)
	}
