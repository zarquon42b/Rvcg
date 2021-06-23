#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/create/resampler.h>
#include <vcg/complex/algorithms/clustering.h>
#include <vcg/simplex/face/distance.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/space/index/grid_static_ptr.h>
#include "typedef.h"
#include <vcg/complex/algorithms/isotropic_remeshing.h>
#include "pointcloud.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
using namespace tri;

using namespace Rcpp;

RcppExport SEXP RisotropicResampling(SEXP vb_, SEXP it_, SEXP TargetLen_, SEXP FeatureAngleDeg_, SEXP MaxSurfDist_, SEXP iterations_, SEXP Adaptive_, SEXP splitFlag_, SEXP collapseFlag_, SEXP swapFlag_, SEXP smoothFlag_,SEXP projectFlag_, SEXP surfDistCheck_) {
  try {

    
    float TargetLen = as<float>(TargetLen_);
    float FeatureAngleDeg = as<float>(FeatureAngleDeg_);
    float MaxSurfDist = as<float>(MaxSurfDist_);
    int iterations = as<int>(iterations_);
    bool Adaptive = as<bool>(Adaptive_);
    bool selectedOnly = false;
    // bool selectedOnly = as<bool>(selectedOnly_);
    bool splitFlag = as<bool>(splitFlag_);
    bool collapseFlag = as<bool>(collapseFlag_);
    bool swapFlag =as<bool>(swapFlag_);
    bool smoothFlag =as<bool>(smoothFlag_);
    bool projectFlag =as<bool>(projectFlag_);
    bool surfDistCheck =as<bool>(surfDistCheck_);
    MyMesh  baseMesh, toProjectCopy;
    int checkit = Rvcg::IOMesh<MyMesh>::RvcgReadR(baseMesh,vb_,it_);
    if (baseMesh.fn==0) {
      ::Rf_error( "This filter requires a mesh with some faces,<br> it does not work on PointSet"); 
      
    }

    tri::Clean<MyMesh>::RemoveDuplicateVertex(baseMesh);
    tri::Clean<MyMesh>::RemoveUnreferencedVertex(baseMesh);
    tri::Allocator<MyMesh>::CompactEveryVector(baseMesh);
    tri::UpdateBounding<MyMesh>::Box(baseMesh);
    baseMesh.face.EnableNormal();
    baseMesh.face.EnableQuality();
    baseMesh.face.EnableFFAdjacency();
    baseMesh.face.EnableVFAdjacency();
    baseMesh.vert.EnableVFAdjacency();
    baseMesh.vert.EnableQuality();
  vcg::tri::Append<MyMesh,MyMesh>::Mesh(toProjectCopy,baseMesh);
  toProjectCopy.face.EnableNormal();
    toProjectCopy.face.EnableQuality();
    toProjectCopy.face.EnableFFAdjacency();
    toProjectCopy.face.EnableVFAdjacency();
    toProjectCopy.vert.EnableVFAdjacency();
    toProjectCopy.vert.EnableQuality();

  //  Point3i volumeDim;
  //Box3f volumeBox = baseMesh.bbox;
  tri::IsotropicRemeshing<MyMesh>::Params params;
		params.SetTargetLen(TargetLen);
		params.SetFeatureAngleDeg(FeatureAngleDeg);


		params.maxSurfDist  = MaxSurfDist;

		params.iter         = iterations;
		params.adapt        = Adaptive;
		params.selectedOnly = selectedOnly;
		params.splitFlag    = splitFlag;
		params.collapseFlag = collapseFlag;
		params.swapFlag     = swapFlag;
		params.smoothFlag   = smoothFlag;
		params.projectFlag  = projectFlag;
		params.surfDistCheck= surfDistCheck;

		
  
		tri::IsotropicRemeshing<MyMesh>::Do(baseMesh, toProjectCopy, params);
		
		return Rvcg::IOMesh<MyMesh>::RvcgToR(baseMesh);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
