#include "typedef.h"
#include "RvcgIO.h" 
#include <Rcpp.h>

using namespace Rcpp;
//#include <wrap/ply/plylib.cpp>
 
RcppExport SEXP Rintersect(SEXP vb_ , SEXP it_, SEXP ioclost_, SEXP normals_, SEXP tol_, SEXP maxtol_, SEXP mindist_)
{
  typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
  ScalarType x,y,z;
  int i;
  VertexIterator vi;
  MyMesh m;
  MyMesh refmesh;
  float tol = Rcpp::as<float>(tol_);
  float maxtol = as<float>(maxtol_);
  Rcpp::NumericMatrix ioclost(ioclost_);
  Rcpp::NumericMatrix normals(normals_);
  int dref =  ioclost.ncol();
  std::vector<float> dis;
  std::vector<float> hitbool;
  // section read from input
  bool mindist = as<bool>(mindist_);
  float t, t1;
  int check = Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
   if (check != 0) {
    Rprintf("%s\n","Warning: mesh has no faces or no vertices, nothing done");
     return Rcpp::List::create(Rcpp::Named("vb") = 0,
			    Rcpp::Named("normals") = 0,
			    Rcpp::Named("hitbool") = 0,
			    Rcpp::Named("dis") = 0
			    );
   }  else {
  //Allocate target
  typedef MyMesh::VertexPointer VertexPointer;
  std::vector<VertexPointer> ivp;
  vcg::tri::Allocator<MyMesh>::AddVertices(refmesh,dref);
  vi=refmesh.vert.begin();
  Point3f normtmp;
  for (i=0; i < dref; i++) {
    x = ioclost(0,i);
    y = ioclost(1,i);
    z = ioclost(2,i);
    (*vi).P() = CoordType(x,y,z);
    x = normals(0,i);
    y = normals(1,i);
    z = normals(2,i);
    (*vi).N() = CoordType(x,y,z);
    ++vi;
  }
  //--------------------------------------------------------------------------------------//
  //
  //                              INITIALIZE SEARCH STRUCTURES
  //
  // Update the FaceProjection flags needed for projection/distance queries
  // Create a static grid (for fast indexing) and fill it 
  //--------------------------------------------------------------------------------------//
  m.face.EnableNormal();
  tri::UpdateBounding<MyMesh>::Box(m);
  tri::UpdateNormal<MyMesh>::PerFaceNormalized(m);//very important !!!
  tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(refmesh);
  float maxDist = m.bbox.Diag();
  float minDist = 1e-10;
  vcg::tri::FaceTmark<MyMesh> mf; 
  mf.SetMesh( &m );
  vcg::RayTriangleIntersectionFunctor<true> FintFunct;
  vcg::face::PointDistanceBaseFunctor<float> PDistFunct;
  TriMeshGrid static_grid;    
  static_grid.Set(m.face.begin(), m.face.end());
  // run search 
  for (i=0; i < refmesh.vn; i++) {
    vcg::Ray3f ray;
    Point3f orig = refmesh.vert[i].P();
    Point3f orig0 = orig;
    Point3f dir = refmesh.vert[i].N();
    Point3f dirOrig = dir;
    Point3f clost = CoordType(0,0,0);
    MyFace* f_ptr; MyFace* f_ptr1;
    t=0; t1=0;
    ray.SetOrigin(orig);
    ray.SetDirection(dir);
    f_ptr = GridDoRay(static_grid, FintFunct, mf, ray, maxDist, t);
    if (! f_ptr) {
      orig0 = orig+CoordType(1e-6,1e-6,1e-6);// add some noise to get eventually hit vertices alon the rays
      ray.SetOrigin(orig0);
      ray.SetDirection(dir);
      f_ptr = GridDoRay(static_grid, FintFunct, mf, ray, maxDist, t);
    }
    if (mindist) {   
      ray.SetDirection(-dir);
      f_ptr1 = GridDoRay(static_grid, FintFunct, mf, ray, maxDist, t1);
      if (! f_ptr1) {
	orig0 = orig+CoordType(1e-6,1e-6,1e-6);// add some noise to get eventually hit vertices alon the rays
	ray.SetOrigin(orig0);
	ray.SetDirection(-dir);
	f_ptr1 = GridDoRay(static_grid, FintFunct, mf, ray, maxDist, t1);
      }
      if ((f_ptr && f_ptr1 && t1 < t) || (!f_ptr && f_ptr1) ) {
	f_ptr = f_ptr1;
	t = -t1;
      } 
    }  
    if (f_ptr && abs(t) < maxtol) {
      if (abs(t) >= tol) {
	clost = refmesh.vert[i].P()+dir*t;//the hit point
	dis.push_back(t);
	hitbool.push_back(1);
      } else {
	dis.push_back(t);
	hitbool.push_back(0);
      }
    } else {
      Point3f& currp = refmesh.vert[i].P();
      f_ptr= GridClosest(static_grid, PDistFunct, mf, currp, maxDist, minDist, clost);
      dis.push_back(minDist);
      hitbool.push_back(0);
    }
  
    int f_i = vcg::tri::Index(m, f_ptr);
    MyMesh::CoordType ti = (m.face[f_i].V(0)->N()+m.face[f_i].V(1)->N()+m.face[f_i].V(2)->N())/3;//the smoothed normal at that point
    ioclost(0,i) = clost[0];
    ioclost(1,i) = clost[1];
    ioclost(2,i) = clost[2];
    normals(0,i) = ti[0];
    normals(1,i) = ti[1];
    normals(2,i) = ti[2];
  }
  return Rcpp::List::create(Rcpp::Named("vb") = ioclost,
			    Rcpp::Named("normals") = normals,
			    Rcpp::Named("hitbool") = hitbool,
			    Rcpp::Named("dis") = dis
			    );
   }
}
    


	
