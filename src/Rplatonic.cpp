
#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
#include<vcg/complex/algorithms/create/platonic.h>

using namespace tri;
using namespace Rcpp;

RcppExport SEXP RSphere(SEXP subdiv_ = wrap(3),SEXP normals_ = wrap(true)) {
  try {
    bool normals = as<bool>(normals_);
    int subdiv = as<int>(subdiv_);
    MyMesh m;
    Sphere(m,subdiv);
    if (normals)
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
    List out = Rvcg::IOMesh<MyMesh>::RvcgToR(m,normals);
    return out;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP RSphericalCap(SEXP angleRad_,SEXP subdiv_ = wrap(3), SEXP normals_ = wrap(true)) {
  try {
    bool normals = as<bool>(normals_);
    int subdiv = as<int>(subdiv_);
    float angleRad = as<float>(angleRad_);
    
    MyMesh m;
    m.vert.EnableVFAdjacency();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();

    SphericalCap(m,angleRad,subdiv);
    if (normals)
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
    List out = Rvcg::IOMesh<MyMesh>::RvcgToR(m,normals);
    return out;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP RTetrahedron(SEXP normals_ = wrap(true)) {
  try {
    bool normals = as<bool>(normals_);
    MyMesh m;
    Tetrahedron(m);
    if (normals)
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
    return Rvcg::IOMesh<MyMesh>::RvcgToR(m,normals);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP RDodecahedron(SEXP normals_ = wrap(true)) {
  bool normals = as<bool>(normals_);
  MyMesh m;
  Dodecahedron(m);
  if (normals)
    tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
  return Rvcg::IOMesh<MyMesh>::RvcgToR(m,normals);
}
RcppExport SEXP ROctahedron(SEXP normals_ = wrap(true)) {
  try {
    bool normals = as<bool>(normals_);
    MyMesh m;
    Octahedron(m);
    if (normals)
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
    return Rvcg::IOMesh<MyMesh>::RvcgToR(m,normals);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP RIcosahedron(SEXP normals_ = wrap(true)) {
  try {
    bool normals = as<bool>(normals_);
    MyMesh m;
    Icosahedron(m);
    if (normals)
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
    return Rvcg::IOMesh<MyMesh>::RvcgToR(m,normals);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}



RcppExport SEXP RHexahedron(SEXP normals_ = wrap(true)) {
  try {
    bool normals = as<bool>(normals_);
    MyMesh m;
    Hexahedron(m);
    if (normals)
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
    return Rvcg::IOMesh<MyMesh>::RvcgToR(m,normals);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP RSquare(SEXP normals_ = wrap(true)) {
  try {
    bool normals = as<bool>(normals_);
    MyMesh m;
    Square(m);
    if (normals)
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
    return Rvcg::IOMesh<MyMesh>::RvcgToR(m,normals);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP RBox(SEXP mesh_,SEXP normals_ = wrap(true)) {
  try {
    List meshraw(mesh_);
    bool normals = as<bool>(normals_);
    MyMesh m, mesh;
    int a = Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(mesh,mesh_);
    //Rvcg::IOMesh<MyMesh>::RvcgReadR(m,meshraw["vb"],meshraw["it"]);
    vcg::Box3<float> bb = mesh.bbox;
    Box(m,bb);
    if (normals)
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
    return Rvcg::IOMesh<MyMesh>::RvcgToR(m,normals);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP RCone(SEXP r1_, SEXP r2_, SEXP h_, SEXP normals_ = wrap(true)) {
  try {
    float r1 = as<float>(r1_);
    float r2 = as<float>(r2_);
     float h = as<float>(h_);
    bool normals = as<bool>(normals_);
    MyMesh m;
    Cone(m,r1,r2,h);
    if (normals)
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
    return Rvcg::IOMesh<MyMesh>::RvcgToR(m,normals);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
