#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/create/resampler.h>
#include <vcg/complex/algorithms/clustering.h>
#include <vcg/simplex/face/distance.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/space/index/grid_static_ptr.h>
#include "typedef.h"
#include "pointcloud.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
using namespace tri;

using namespace Rcpp;

RcppExport SEXP RuniformResampling(SEXP vb_, SEXP it_, SEXP voxelSize_, SEXP offsetThr_, SEXP discretizeFlag_, SEXP multiSampleFlag_, SEXP absDistFlag_, SEXP mergeCloseVert_, SEXP silent_) {
  try {
  float voxelSize = as<float>(voxelSize_);
  float offsetThr = as<float>(offsetThr_);
  bool discretizeFlag = as<bool>(discretizeFlag_);
  bool multiSampleFlag = as<bool>(multiSampleFlag_);
  bool absDistFlag = as<bool>(absDistFlag_);
  bool mergeCloseVert =as<bool>(mergeCloseVert_);
  bool silent =as<bool>(silent_);
  MyMesh m, baseMesh, offsetMesh;
  int checkit = Rvcg::IOMesh<MyMesh>::RvcgReadR(baseMesh,vb_,it_);
  if (baseMesh.fn==0) {
    ::Rf_error( "This filter requires a mesh with some faces,<br> it does not work on PointSet"); 

  }
  tri::UpdateBounding<MyMesh>::Box(baseMesh);
  baseMesh.face.EnableNormal();
  Point3i volumeDim;
  Box3f volumeBox = baseMesh.bbox;
  volumeBox.Offset(volumeBox.Diag()/10.0f+offsetThr);
  
  BestDim(volumeBox , voxelSize, volumeDim );
  if (!silent) {
    Rprintf("     Resampling mesh using a volume of %i x %i x %i\n",volumeDim[0],volumeDim[1],volumeDim[2]);
    Rprintf("     VoxelSize is %f, offset is %f\n", voxelSize,offsetThr);
    Rprintf("     Mesh Box is %f %f %f\n",baseMesh.bbox.DimX(),baseMesh.bbox.DimY(),baseMesh.bbox.DimZ() );
  }
  tri::Resampler<MyMesh,MyMesh>::Resample(baseMesh, offsetMesh, volumeBox, volumeDim, voxelSize*3.5, offsetThr,discretizeFlag,multiSampleFlag,absDistFlag);
  if (mergeCloseVert) {
    float mergeThr =offsetMesh.bbox.Diag()/10000.0f;
    int total = tri::Clean<MyMesh>::MergeCloseVertex(offsetMesh,mergeThr);
    if (!silent)
      Rprintf("\nSuccessfully merged %d vertices with a distance lower than %f\n", total, mergeThr);
  }
  vcg::tri::Allocator< MyMesh >::CompactVertexVector(offsetMesh);
  vcg::tri::Allocator< MyMesh >::CompactFaceVector(offsetMesh);
  tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(offsetMesh);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(offsetMesh);
  Rcpp::NumericMatrix vbout(3,offsetMesh.vn), normals(3,offsetMesh.vn);
  Rcpp::IntegerMatrix itout(3,offsetMesh.fn);
  SimpleTempData<MyMesh::VertContainer,int> indiceout(offsetMesh.vert);
  VertexIterator vi;
  FaceIterator fi;
 vi=offsetMesh.vert.begin();
  for (int i=0;  i < offsetMesh.vn; i++) {
    indiceout[vi] = i;
    vbout(0,i) = (*vi).P()[0];
    vbout(1,i) = (*vi).P()[1];
    vbout(2,i) = (*vi).P()[2];
    normals(0,i) = (*vi).N()[0];
    normals(1,i) = (*vi).N()[1];
    normals(2,i) = (*vi).N()[2];
    ++vi;
  }
  FacePointer fp;
  
  fi=offsetMesh.face.begin();
  int j = 0;
  for (int i=0; i < offsetMesh.fn; i++) {
    fp=&(*fi);
    itout(0,i) = indiceout[fp->cV(0)]+1;
    itout(1,i) = indiceout[fp->cV(1)]+1;
    itout(2,i) = indiceout[fp->cV(2)]+1;
    ++fi;
  }
  
  return Rcpp::List::create(Rcpp::Named("vb") = vbout,
			    Rcpp::Named("it") = itout,
			    Rcpp::Named("normals")=normals);
  
  return wrap(0);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
