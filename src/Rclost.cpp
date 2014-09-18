#include "typedef.h"
#include "pointcloud.h"
#include "RvcgIO.h"
#include <Rcpp.h>
#include <cmath>

using namespace tri;

using namespace Rcpp;

RcppExport SEXP Rclost(SEXP vb_ , SEXP it_, SEXP ioclost_, SEXP sign_, SEXP borderchk_, SEXP barycentric_, SEXP smooth_)
{
  try {
    typedef vcg::SpatialHashTable<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid; 
    //typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
    bool signo = as<bool>(sign_);
    bool borderchk = as<bool>(borderchk_);
    bool barycentric = as<bool>(barycentric_);
    bool smooth = as<bool>(smooth_);
    Rcpp::NumericMatrix ioclost(ioclost_);
    int i;
    MyMesh m;
    PcMesh refmesh;
    PcMesh outmesh;
    MyMesh::CoordType baryco;
    // section read from input
    int checkit = Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    if (checkit == 1) {
      ::Rf_error("target mesh has no faces, nothing done");
      
    } else if (checkit >= 0) {
      Rvcg::IOMesh<PcMesh>::RvcgReadR(refmesh, ioclost_); 
      m.face.EnableNormal();
 
      tri::UpdateBounding<MyMesh>::Box(m);
      tri::UpdateNormal<MyMesh>::PerFaceNormalized(m);//very important !!!
      //tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
      if (smooth) {
	tri::Smooth<MyMesh>::VertexNormalLaplacian(m,2,false);
	tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
      }
      float maxDist = m.bbox.Diag()*2;
      float minDist = 1e-10;
      vcg::tri::FaceTmark<MyMesh> mf; 
      mf.SetMesh( &m );
      vcg::face::PointDistanceBaseFunctor<float> PDistFunct;
      TriMeshGrid static_grid;    
      static_grid.Set(m.face.begin(), m.face.end());
      if (borderchk) { //update Border flags
	m.vert.EnableVFAdjacency();
	m.face.EnableFFAdjacency();
	m.face.EnableVFAdjacency();
	tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
	tri::UpdateSelection<MyMesh>::FaceFromBorderFlag(m);
      }
      //setup return structure
      Rcpp::NumericMatrix normals(3,refmesh.vn), barycoord(3,refmesh.vn);
      Rcpp::IntegerVector border(refmesh.vn), faceptr(refmesh.vn);
      Rcpp::NumericVector dis(refmesh.vn);
      //index faces
      SimpleTempData<MyMesh::FaceContainer,int> indices(m.face);
      FaceIterator fi=m.face.begin();
      for (i=0; i < m.fn; i++) {
	indices[fi] = i;
	++fi;
      }
      vcg::tri::Append<PcMesh,PcMesh>::Mesh(outmesh,refmesh);
      PcMesh::CoordType tt;
      for(i=0; i < refmesh.vn; i++) {
	border[i] = 0;
	Point3f& currp = refmesh.vert[i].P();
	Point3f& clost = outmesh.vert[i].P();
	MyFace* f_ptr= GridClosest(static_grid, PDistFunct, mf, currp, maxDist, minDist, clost);
	if (f_ptr) {
	  if (borderchk) {
	    if ((*f_ptr).IsS())
	      border[i] = 1;
	  }
	  faceptr[i] = indices[f_ptr];
	  int f_i = vcg::tri::Index(m, f_ptr);
	  tt = currp*0;
	
	  for (int j=0; j <3;j++) {
	    if (&(m.face[f_i].V(j)->N())) {
	      Point3f vdist = m.face[f_i].V(j)->P() - clost;
	      float weight = sqrt(vdist.dot(vdist));
	      if (weight > 0)
		weight = 1/weight;
	      else 
		weight = 1e12;
	    
	      tt +=(m.face[f_i].V(j)->N()*weight);
	    }
	  }
	  if (barycentric) {
	    baryco = currp*0;
	    InterpolationParameters<MyFace,ScalarType>(*f_ptr,f_ptr->N(),clost,baryco);
	  }
	}
	float vl = sqrt(tt.dot(tt));
	if (vl > 0) {//check for zero length normals
	  tt=tt/vl;
	}   
    
	dis[i] = minDist;
	if (signo) {
	  Point3f dif = clost - currp;
	  float sign = dif.dot(tt);	
	  if (sign < 0)
	    dis[i] = -dis[i] ;
	}
	//write back output
	ioclost(0,i) =clost[0];
	ioclost(1,i) =clost[1];
	ioclost(2,i) =clost[2];
	normals(0,i) = tt[0];
	normals(1,i) = tt[1];    
	normals(2,i) = tt[2];
	if(barycentric) {
	  barycoord(0,i) = baryco[0];
	  barycoord(1,i) = baryco[1];
	  barycoord(2,i) = baryco[2];
	}
      
      }
      return Rcpp::List::create(Rcpp::Named("ioclost") = ioclost,
				Rcpp::Named("barycoord") = barycoord,
				Rcpp::Named("normals") = normals,
				Rcpp::Named("border") = border, 
				Rcpp::Named("distance") = dis,
				Rcpp::Named("faceptr") = faceptr

				);
    } else {
        
      return wrap(1);
    }
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

