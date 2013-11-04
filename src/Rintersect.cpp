#include "typedef.h"
#include "RvcgIO.h" 
#include <Rcpp.h>

using namespace Rcpp;
//#include <wrap/ply/plylib.cpp>
 
RcppExport SEXP Rintersect(SEXP _vb , SEXP _it, SEXP _ioclost, SEXP _normals, SEXP _tol)
{
  typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
  ScalarType x,y,z;
  int i;
  VertexIterator vi;
  MyMesh m;
  MyMesh refmesh;
  MyMesh outmesh;
  float tol = Rcpp::as<float>(_tol);	
  Rcpp::NumericMatrix ioclost(_ioclost);
  Rcpp::NumericMatrix normals(_normals);
  int dref =  ioclost.ncol();
  std::vector<float> dis;
  std::vector<float> hitbool;
  // section read from input
   
  float t;
  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,_vb,_it);
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
  vcg::tri::Append<MyMesh,MyMesh>::Mesh(outmesh,refmesh);
  tri::UpdateBounding<MyMesh>::Box(m);
  tri::UpdateNormal<MyMesh>::PerFaceNormalized(m);//very important !!!
  tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(refmesh);
  float maxDist = m.bbox.Diag();
  float minDist = 1e-10;
  vcg::tri::FaceTmark<MyMesh> mf; 
  vcg::tri::VertTmark<MyMesh> mv;
  mf.SetMesh( &m );
  mv.SetMesh( &m );
  vcg::RayTriangleIntersectionFunctor<true> FintFunct;
  vcg::face::PointDistanceBaseFunctor<float> PDistFunct;
  vcg::vertex::PointNormalDistanceFunctor<MyVertex> VDistFunct;
  TriMeshGrid static_grid;    
  static_grid.Set(m.face.begin(), m.face.end());
  // run search 
  for (i=0; i < refmesh.vn; i++) {
    vcg::Ray3f ray;
    Point3f orig = refmesh.vert[i].P();
    Point3f dir = refmesh.vert[i].N();
    Point3f dirOrig = dir;
    Point3f clost = CoordType(0,0,0);
    MyFace* f_ptr;
    t=0;
    ray.SetOrigin(orig);
    ray.SetDirection(dir);
    f_ptr = GridDoRay(static_grid, FintFunct, mf, ray, maxDist, t);
      
    if (f_ptr) {
      if (t >= tol) {
	clost = refmesh.vert[i].P()+dir*t;//the hit point
	dis.push_back(t);
	hitbool.push_back(1);
      } else {
	dis.push_back(t);
	hitbool.push_back(0);
      }
    } else {
      orig = orig+CoordType(1e-6,1e-6,1e-6);// add some noise to get eventually hit vertices alon the rays
      ray.SetOrigin(orig);
      ray.SetDirection(dir);
      f_ptr = GridDoRay(static_grid, FintFunct, mf, ray, maxDist, t);
      if (f_ptr) {
	if (t >= tol) {
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
    


	
