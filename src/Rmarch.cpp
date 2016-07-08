#include "typedef.h"
#include "Volume.h"

//#include "RvcgIO.h"
#include <RcppArmadillo.h>
//#include <RcppArmadillo.h>


//using namespace Rcpp;
//using namespace std;
//using namespace vcg;
using namespace tri;
using namespace Rcpp;
//using namespace arma;
#include <vcg/complex/complex.h>
#include <vcg/math/perlin_noise.h>
#include <vcg/complex/algorithms/create/marching_cubes.h>
//#include <vcg/complex/algorithms/create/extended_marching_cubes.h>
#include <vcg/complex/algorithms/create/mc_trivial_walker.h>
#include <wrap/io_trimesh/export_ply.h>
//#include "simple_volume.h"
//using namespace std;
//#include "Voxel.h"

RcppExport SEXP RMarchC(SEXP array_, SEXP thresh_) {
  try {
    IntegerVector tmparr(array_);
    IntegerVector arrayDims = tmparr.attr("dim");
   
    std::vector<float> vecArray = as<std::vector<float> >(array_);
    
    //bin/IntegerVector vecArray(array_);
    double thresh= as<double>(thresh_);
    
  
MyMesh m;
VertexIterator vi;
FaceIterator fi;
// allocate mesh and fill it
/*m.vert.EnableVFAdjacency();
  m.face.EnableFFAdjacency();
  m.face.EnableVFAdjacency();*/
int i,j,k;

typedef MySimpleVolume<MySimpleVoxel> MyVolume;
MyVolume	volume;
typedef vcg::tri::TrivialWalker<MyMesh,MyVolume>	MyWalker;
typedef vcg::tri::MarchingCubes<MyMesh, MyWalker>	MyMarchingCubes;
MyWalker walker;
volume.Init(Point3i(arrayDims[0],arrayDims[1],arrayDims[2]));
 for(i=0;i < arrayDims[0];i++) {
   for(j=0;j<arrayDims[1];j++) {
     for(k=0;k< arrayDims[2];k++) {
       int tmpval = vecArray[i+j*arrayDims[0]+k*(arrayDims[0]*arrayDims[1])];
       /*if (tmpval >= lower && tmpval <= upper)
	 volume.Val(i,j,k)=tmpval;
	 else*/ 
	 volume.Val(i,j,k)=tmpval;
     }
   }
 }
  //write back
/*volume.Init(Point3i(64,64,64));
  for(int i=0;i<64;i++)
    for(int j=0;j<64;j++)
      for(int k=0;k<64;k++)
      volume.Val(i,j,k)=(j-32)*(j-32)+(k-32)*(k-32)  + i*10*(float)math::Perlin::Noise(i*.2,j*.2,k*.2);*/
MyMarchingCubes	mc(m, walker);
walker.BuildMesh<MyMarchingCubes>(m, volume, mc, thresh);
  vcg::tri::Allocator< MyMesh >::CompactVertexVector(m);
  vcg::tri::Allocator< MyMesh >::CompactFaceVector(m);
  tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
  SimpleTempData<MyMesh::VertContainer,int> indiceout(m.vert);
  Rcpp::NumericMatrix vbout(3,m.vn), normals(3,m.vn);
  Rcpp::IntegerMatrix itout(3,m.fn);
  
  vi=m.vert.begin();
  for (i=0;  i < m.vn; i++) {
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
  
  fi=m.face.begin();
  j = 0;
  for (i=0; i < m.fn; i++) {
    fp=&(*fi);
    itout(0,i) = indiceout[fp->cV(0)]+1;
    itout(1,i) = indiceout[fp->cV(1)]+1;
    itout(2,i) = indiceout[fp->cV(2)]+1;
    ++fi;
  }
  //delete &walker;
  //delete &volume;
  return Rcpp::List::create(Rcpp::Named("vb") = vbout,
			    Rcpp::Named("it") = itout,
			    Rcpp::Named("normals") = normals
			      );
 } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }

}
 

    
