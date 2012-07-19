// Author: Stefan Schlager
// Date: 15 September 2010
/*
#include <string.h>
#include <vector>
using namespace std;
#include <stdio.h>
#include <cstddef>

// VCG headers for triangular mesh processing
#include<vcg/simplex/edge/base.h>
#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/face/base.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/edges.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/intersection.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/smooth.h>
#include<vcg/complex/allocate.h>
#include <wrap/callback.h>
#include <vcg/complex/append.h>

// VCG File Format Importer/Exporter
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/update/color.h>*/
#include <../typedef.h>
#include <R.h>
#include <Rdefines.h> 
//#include <Rcpp.h>

//using namespace std;
 extern "C" {

   void testfield(double *data ,int *dim)
  {
    
    ScalarType x,y,z;
    int i;
    const int d = *dim;
    MyMesh m;
    int n = 5;
    vcg::tri::Allocator<MyMesh>::AddVertices(m,d);
    
    VertexPointer ivp[n];
    VertexIterator vi=m.vert.begin();
    for (i=0; i<d; i++) 
    {
      
      x = data[i*3];
      y = data[i*3+1];
      z=  data[i*3+2];
      (*vi).P() = CoordType(x,y,z);
      ++vi;
    }
    
    // SEXP tmp;
    //(tmp = allocVector(REALSXP, 27));
 
    //SET_LENGTH(tmp,27);
    //PROTECT(tmp = NEW_NUMERIC(27)) ; 
    // update output
       double tmp[27];
     vi=m.vert.begin();
    for (i = 0;i<8;i++)
      {	
	
	data[i*3] = (*vi).P()[0];
	data[i*3+1];
	data[i*3+2]+1;	
	++vi;
	}
    for (i = 9;i < d;i++)
      {data[i*3]=0;
	  data[i*3+1]=0;
	  data[i*3+2]=0;
      }
    //delete data;
    //data = tmp;
    //SEXP tt = 5;
    // Rcpp::wrap(tmp);
    
  }


  void Rsmooth(double *vb ,int *dim, int *it, int *dimit, int *iteration, int *stype, double *normals)
  {
    typedef MyMesh::CoordType CoordType;
    typedef MyMesh::ScalarType ScalarType;
    ScalarType x,y,z;
    int i;
    const int d = *dim;
    const int faced = *dimit;
    int iter = *iteration;
    int method = *stype;
    MyMesh m;
    //int n = 5;
    vcg::tri::Allocator<MyMesh>::AddVertices(m,d);
    vcg::tri::Allocator<MyMesh>::AddFaces(m,faced);
    VertexPointer ivp[d];
    VertexIterator vi=m.vert.begin();
    for (i=0; i<d; i++) 
      {
	ivp[i]=&*vi;
	x = vb[i*3];
	y = vb[i*3+1];
	z=  vb[i*3+2];
	(*vi).P() = CoordType(x,y,z);
	++vi;
      }
    int itx,ity,itz;
    FaceIterator fi=m.face.begin();
    for (i=0; i < faced; i++) 
      {
	itx = it[i*3];
	ity = it[i*3+1];
	itz = it[i*3+2];
	(*fi).V(0)=ivp[itx];
	(*fi).V(1)=ivp[ity];
	(*fi).V(2)=ivp[itz];
	++fi;
      }
    if (method == 0)
      {
	tri::Smooth<MyMesh>::VertexCoordTaubin(m,iter,0.5,-0.53);
      }
    else if (method == 1)
      {
	tri::Smooth<MyMesh>::VertexCoordLaplacian(m,iter);
      }
    else if (method == 2)
      {
	tri::Smooth<MyMesh>::VertexCoordLaplacianHC(m,iter);
      }
    tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(m);
    tri::UpdateNormals<MyMesh>::NormalizeVertex(m);
   
    //write back output
    vi=m.vert.begin();
    for (i=0; i<d; i++) 
      {
	vb[i*3] = (*vi).P()[0];
	vb[i*3+1] = (*vi).P()[1];
	vb[i*3+2] = (*vi).P()[2];
	normals[i*3] = (*vi).N()[0];
	normals[i*3+1] = (*vi).N()[1];
	normals[i*3+2] = (*vi).N()[2];
	++vi;
      }
  }
  
   

}
