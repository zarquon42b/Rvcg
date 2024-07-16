// Author: Stefan Schlager
// Date: 15 September 2010

#include <vcg/complex/algorithms/implicit_smooth.h>
#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>

using namespace Rcpp;

RcppExport SEXP Rsmooth(SEXP vb_, SEXP it_, SEXP iteration_, SEXP method_, SEXP lambda_,  SEXP mu_, SEXP delta_)
{
  try {
    int i;
    MyMesh m;
    VertexIterator vi;
    FaceIterator fi;
    //set up parameters
    int iter = Rcpp::as<int>(iteration_);
    int method = Rcpp::as<int>(method_);
    float lambda =  Rcpp::as<float>(lambda_);
    float mu =  Rcpp::as<float>(mu_);
    ScalarType delta = Rcpp::as<double>(delta_);
    //allocate mesh and fill it
    Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);

    if (method == 0) {
      tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
      unsigned int cnt=tri::UpdateSelection<MyMesh>::VertexFromFaceStrict(m);
      tri::Smooth<MyMesh>::VertexCoordTaubin(m, iter, lambda, mu, cnt>0);
    } else if (method == 1) {
      tri::Smooth<MyMesh>::VertexCoordLaplacian(m, iter);
    } else if (method == 2) {
      tri::UpdateSelection<MyMesh>::FaceAll(m);
      tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
      unsigned int cnt=tri::UpdateSelection<MyMesh>::VertexFromFaceStrict(m);
      tri::Smooth<MyMesh>::VertexCoordLaplacianHC(m, iter,cnt>0);
    } else if (method == 3) {
      tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
      tri::UpdateFlags<MyMesh>::FaceClearB(m);
      tri::Smooth<MyMesh>::VertexCoordScaleDependentLaplacian_Fujiwara(m,iter,delta);
    } else if (method == 4) {
      tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
      tri::UpdateFlags<MyMesh>::FaceClearB(m);
      tri::Smooth<MyMesh>::VertexCoordLaplacianAngleWeighted(m,iter,delta);
    }
    else if (method == 5) {
      tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
      tri::UpdateFlags<MyMesh>::FaceClearB(m);
      tri::Smooth<MyMesh>::VertexCoordPlanarLaplacian(m, iter, delta);
    }
    vcg::tri::Allocator<MyMesh>::CompactVertexVector(m);
    vcg::tri::Allocator<MyMesh>::CompactFaceVector(m);
    tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
    tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
    Rcpp::NumericMatrix vb(3, m.vn);
    Rcpp::NumericMatrix normals(3, m.vn);
    Rcpp::IntegerMatrix itout(3, m.fn);
    //write back output
    SimpleTempData<MyMesh::VertContainer,int>indices(m.vert);

    // write back updated mesh
    vi=m.vert.begin();
    for (i=0; i < m.vn; i++) {
      indices[vi] = i;
      if( ! vi->IsD() ) {
	vb(0,i) = (*vi).P()[0];
	vb(1,i) = (*vi).P()[1];
	vb(2,i) = (*vi).P()[2];
	normals(0,i) = (*vi).N()[0];
	normals(1,i) = (*vi).N()[1];
	normals(2,i) = (*vi).N()[2];
      }
      ++vi;
    }

    FacePointer fp;
    fi=m.face.begin();
    for (i=0; i < m.fn; i++) {
      fp=&(*fi);
      if( ! fp->IsD() ) {
	itout(0,i) = indices[fp->cV(0)]+1;
	itout(1,i) = indices[fp->cV(1)]+1;
	itout(2,i) = indices[fp->cV(2)]+1;
      }
      ++fi;
    }
    return Rcpp::List::create(Rcpp::Named("vb") = vb,
			      Rcpp::Named("normals") = normals,
			      Rcpp::Named("it") = itout
			      );

  } catch (std::exception& e) {
    forward_exception_to_r( e );
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP RsmoothImplicit(
    SEXP vb_, SEXP it_, SEXP lambda_, SEXP useMassMatrix_, SEXP fixBorder_,
    SEXP useCotWeight_, SEXP degree_, SEXP lapWeight_, SEXP SmoothQ_)
{
  try {
    int i;
    MyMesh m;
    VertexIterator vi;
    FaceIterator fi;
    //set up parameters
    ScalarType lambda = Rcpp::as<float>(lambda_);
    //the use of mass matrix to keep the mesh close to its original position
    //(weighted per area distributed on vertices)
    bool useMassMatrix = Rcpp::as<bool>(useMassMatrix_);
    //this bool is used to fix the border vertices of the mesh or not
    bool fixBorder = Rcpp::as<bool>(fixBorder_);
    //this bool is used to set if cotangent weight is used, this flag to false means uniform laplacian
    bool useCotWeight = Rcpp::as<bool>(useCotWeight_);
    //the set of fixed vertices
    // std::vector<int> FixedV;
    //the set of faces for barycentric constraints
    // std::vector<FaceConstraint> ConstrainedF;
    //the degree of laplacian
    int degree = Rcpp::as<int>(degree_);
    //use this weight for the laplacian when the cotangent one is not used
    ScalarType lapWeight = Rcpp::as<float>(lapWeight_);
    //this is to say if we smooth the positions or the quality
    bool SmoothQ = Rcpp::as<bool>(SmoothQ_);


    //allocate mesh and fill it
    Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);


    vcg::ImplicitSmoother<MyMesh>::Parameter par;
    par.lambda = lambda;
    par.useMassMatrix = useMassMatrix;
    par.fixBorder = fixBorder;
    par.useCotWeight = useCotWeight;
    par.degree = degree;
    par.lapWeight = lapWeight;
    par.SmoothQ = SmoothQ;

    vcg::ImplicitSmoother<MyMesh>::Compute(m, par);

    Rcpp::checkUserInterrupt();

    vcg::tri::Allocator<MyMesh>::CompactVertexVector(m);
    vcg::tri::Allocator<MyMesh>::CompactFaceVector(m);
    tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
    tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
    Rcpp::NumericMatrix vb(3, m.vn);
    Rcpp::NumericMatrix normals(3, m.vn);
    Rcpp::IntegerMatrix itout(3, m.fn);
    //write back output
    SimpleTempData<MyMesh::VertContainer,int>indices(m.vert);

    Rcpp::checkUserInterrupt();
    // write back updated mesh
    vi=m.vert.begin();
    for (i=0; i < m.vn; i++) {
      indices[vi] = i;
      if( ! vi->IsD() ) {
        vb(0,i) = (*vi).P()[0];
        vb(1,i) = (*vi).P()[1];
        vb(2,i) = (*vi).P()[2];
        normals(0,i) = (*vi).N()[0];
        normals(1,i) = (*vi).N()[1];
        normals(2,i) = (*vi).N()[2];
      }
      ++vi;
    }

    FacePointer fp;
    fi=m.face.begin();
    for (i=0; i < m.fn; i++) {
      fp=&(*fi);
      if( ! fp->IsD() ) {
        itout(0,i) = indices[fp->cV(0)]+1;
        itout(1,i) = indices[fp->cV(1)]+1;
        itout(2,i) = indices[fp->cV(2)]+1;
      }
      ++fi;
    }
    return Rcpp::List::create(Rcpp::Named("vb") = vb,
                              Rcpp::Named("normals") = normals,
                              Rcpp::Named("it") = itout
    );

  } catch (std::exception& e) {
    forward_exception_to_r( e );
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

