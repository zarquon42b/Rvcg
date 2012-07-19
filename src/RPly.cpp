/*#include <string.h>
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
#include <wrap/ply/plylib.cpp>
#include <vcg/container/simple_temporary_data.h>
#include <wrap/io_trimesh/import.h>
#include <string.h>
  
  
extern "C" {

  void RPlyRead(char **filename, double *vb ,int *dim, int *it, int *dimit, double *normals, int *getNorm, int *updNorm, double *quality)
  {

    ScalarType x,y,z;
    int i;
    MyMesh m;
       // section read from input
    int d = *dim;
    int faced = *dimit;
    //char file = **filename;
    char file[256];
    strcpy(file, *filename);
    
    int importNorm = *getNorm;
    int updateNorm = *updNorm;
    //load file
    int err2 = tri::io::ImporterPLY<MyMesh>::Open(m,file);
    if(err2) {
     printf("Error in reading %s: '%s'\n",file,tri::io::Importer<MyMesh>::ErrorMsg(err2));
     //exit(-1);  
     }
    //printf("%i",err2);
    if (err2 == 0)
      {
	if (updateNorm == 1)
	  {
	    tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(m);
	    tri::UpdateNormals<MyMesh>::NormalizeVertex(m);
	  }
    //--------------------------------------------------------------------------------------//
    //
    //                                   WRITE BACK
    // Create meshes,
    // Update the bounding box and initialize max search distance
    // Remove duplicates and update mesh properties
    //--------------------------------------------------------------------------------------//
	SimpleTempData<typename MyMesh::VertContainer,int> indices(m.vert);
	
	//VertexPointer ivp[d];
	if (m.vn > 0)
	  {
	    VertexIterator vi=m.vert.begin();
	    
	    for (i=0;  i < m.vn; i++) 
	      {
		indices[vi] = i;//important: updates vertex indices
		//	ivp[i]=&*vi;
		vb[i*3] = (*vi).P()[0];
		vb[i*3+1] = (*vi).P()[1];
		vb[i*3+2] = (*vi).P()[2];
		if (importNorm == 1)
		  {
		    normals[i*3] = (*vi).N()[0];
		    normals[i*3+1] = (*vi).N()[1];
		    normals[i*3+2] = (*vi).N()[2];
		  }
		++vi;
	      }
	  }
	
	FacePointer fp;
	int vv[3];
	*dim = m.vn;
	FaceIterator fi=m.face.begin();
	faced=m.fn;
	if (m.fn > 0)
	  {
	    for (i=0; i < faced;i++) 
	      {
		fp=&(*fi);
		if( ! fp->IsD() )
		  {
		    vv[0]=indices[fp->cV(0)];
		    vv[1]=indices[fp->cV(1)];
		    vv[2]=indices[fp->cV(2)];
		    it[i*3]=vv[0];
		    it[i*3+1]=vv[1];
		    it[i*3+2]=vv[2];
		    ++fi;
		  }
	      }
	  }
	*dimit=m.fn;
      }
  }
}
