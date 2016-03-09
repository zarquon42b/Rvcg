#include "typedef.h"
#include "pointcloud.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
#include <cmath>

using namespace tri;

using namespace Rcpp;

RcppExport SEXP Rclost(SEXP mesh_, SEXP ioclost_, SEXP sign_, SEXP borderchk_, SEXP barycentric_, SEXP smooth_,SEXP tol_, SEXP facenormals_ = wrap(true)) {
  try {
    typedef vcg::SpatialHashTable<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid; 
    //typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
    bool signo = as<bool>(sign_);
    bool borderchk = as<bool>(borderchk_);
    bool barycentric = as<bool>(barycentric_);
    bool smooth = as<bool>(smooth_);
    bool facenormals = as<bool>(facenormals_);
    float tol = as<float>(tol_);
    
    int i;
    MyMesh m;
    PcMesh refmesh;
    PcMesh outmesh;
    MyMesh::CoordType baryco;
    // section read from input
    int checkit = Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(m,mesh_);
    if (checkit == 1)
      ::Rf_error("target mesh has no faces, nothing done");
    
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
    if (tol > 0) 
      maxDist = tol;
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
    
    //setup output datastructures
    arma::mat ioclost(4,refmesh.vn), normals(4,refmesh.vn);
    ioclost.fill(1);normals.fill(1);
    arma::mat barycoord;
    arma::ivec border, faceptr(refmesh.vn);
    arma::vec distances(refmesh.vn);
    if (barycentric)
      barycoord.resize(3,refmesh.vn);
    if (borderchk) {
      border.resize(refmesh.vn);
      border.fill(0);
    }
    //index faces
    SimpleTempData<MyMesh::FaceContainer,int> indices(m.face);
    FaceIterator fi=m.face.begin();
    for (i=0; i < m.fn; i++) {
      indices[fi] = i;
      ++fi;
    }
    vcg::tri::Append<PcMesh,PcMesh>::Mesh(outmesh,refmesh);
    PcMesh::CoordType vertexnormal;
    for(i=0; i < refmesh.vn; i++) {
      
      Point3f& currp = refmesh.vert[i].P();
      Point3f& clost = outmesh.vert[i].P();
      MyFace* f_ptr= GridClosest(static_grid, PDistFunct, mf, currp, maxDist, minDist, clost);
      if (f_ptr) {
	if (borderchk) {
	  if ((*f_ptr).IsS())
	    border[i] = 1;
	}
	faceptr[i] = indices[f_ptr]+1;
	int f_i = vcg::tri::Index(m, f_ptr);
	vertexnormal = currp*0;
	if (facenormals) {
	  vertexnormal = m.face[f_i].N();
	} else {
	  for (int j=0; j <3;j++) {
	    //if (&(m.face[f_i].V(j)->N())) {
	    Point3f vdist = m.face[f_i].V(j)->P() - clost;
	    float weight = sqrt(vdist.dot(vdist));
	    if (weight > 0)
	      weight = 1/weight;
	    else 
	      weight = 1e12;
	    vertexnormal +=(m.face[f_i].V(j)->N()*weight);
	  }
	}
	if (barycentric) {
	  baryco = currp*0;
	  InterpolationParameters<MyFace,ScalarType>(*f_ptr,f_ptr->N(),clost,baryco);
	}
	float vl = sqrt(vertexnormal.dot(vertexnormal));
	if (vl > 0) {//check for zero length normals
	  vertexnormal=vertexnormal/vl;
	  distances[i] = minDist;
	  if (signo) {
	    Point3f dif = clost - currp;
	    float sign = dif.dot(vertexnormal);	
	    if (sign < 0)
	      distances[i] = -distances[i] ;
	  }
	}
      } else {
	double mynan = std::nan("1");
	distances[i] = mynan;
      }
      //write back output
      for (int j = 0; j < 3; j++) {
	ioclost(j,i) =clost[j];
	normals(j,i) = vertexnormal[j];
	if(barycentric) 
	  barycoord(j,i) = baryco[j];
      }
    }
      List out = Rcpp::List::create(Rcpp::Named("vb") = ioclost,
				    Rcpp::Named("it")=wrap(1),
				    Rcpp::Named("normals") = normals,
				    Rcpp::Named("quality") = NumericVector(distances.begin(),distances.end()),
				    Rcpp::Named("faceptr") = NumericVector(faceptr.begin(),faceptr.end())
				    );
      if (barycentric)
	out["barycoords"] = barycoord;
      if (borderchk)
	out["border"] = NumericVector(border.begin(),border.end());

      out.attr("class") = "mesh3d";
      return out;
     
    
    } catch (std::exception& e) {
      ::Rf_error( e.what());
    } catch (...) {
      ::Rf_error("unknown exception");
    }
  }

