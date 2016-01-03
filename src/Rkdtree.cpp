#include "typedef.h"
#include "pointcloud.h"
#include "RvcgIO.h"
#include <Rcpp.h>
#include <RvcgKD.h>
#include <Rconfig.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace tri;
using namespace Rcpp;

RcppExport SEXP Rkdtree(SEXP vb0_, SEXP vb1_, SEXP k_ ,SEXP nofP_= wrap(16),SEXP mDepth_= wrap(64),SEXP threads_=wrap(1)) {
  try {
    int k = as<int>(k_);
    int threads = as<int>(threads_);
    unsigned int nofP = as<unsigned int >(nofP_);
    unsigned int mDepth = as<unsigned int >(mDepth_);
    typedef pair<float,int> mypair;
    PcMesh target, query;
    Rvcg::IOMesh<PcMesh>::RvcgReadR(target, vb0_);  
    Rvcg::IOMesh<PcMesh>::RvcgReadR(query, vb1_);
 
    List out = Rvcg::KDtree< PcMesh, PcMesh >::KDtreeIO(target, query, k,nofP, mDepth,threads);
    return out;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
  
}

RcppExport SEXP RclosestKD(SEXP vb_, SEXP it_, SEXP ioclost_, SEXP itclost_, SEXP k_, SEXP sign_, SEXP smooth_, SEXP barycentric_, SEXP borderchk_, SEXP nofP_= wrap(16),SEXP mDepth_= wrap(64),SEXP angdev_=wrap(0), SEXP wnorm_=wrap(true), SEXP facenormals_ = wrap(true), SEXP threads_=wrap(1)) {
  try {
    bool smooth = as<bool>(smooth_);
    bool barycentric = as<bool>(barycentric_);
    bool borderchk = as<bool>(borderchk_);
    bool wnorm = as<bool>(wnorm_);
    bool facenormals = as<bool>(facenormals_);
    unsigned int nofP = as<unsigned int >(nofP_);
    unsigned int mDepth = as<unsigned int >(mDepth_);
    int threads = as<int>(threads_);
    int k = as<int>(k_);
    bool sign = as<bool>(sign_);
    MyMesh target;
    PcMesh  bary;
    MyMesh query;
    int checkit = Rvcg::IOMesh<MyMesh>::RvcgReadR(target,vb_,it_);
    double angdev = as<double>(angdev_);
    target.face.EnableNormal();
    checkit = Rvcg::IOMesh<MyMesh>::RvcgReadR(query, ioclost_,itclost_);
    if (angdev > 0) {
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(query);
    }
    tri::UpdateNormal<MyMesh>::PerFaceNormalized(target);
    tri::UpdateNormal<MyMesh>::PerVertexNormalized(target);
    if (smooth) {
      tri::Smooth<MyMesh>::VertexNormalLaplacian(target,2,false);
      tri::UpdateNormal<MyMesh>::NormalizePerVertex(target);
    }
    if (borderchk) { //update Border flags
      tri::UpdateFlags<MyMesh>::FaceBorderFromNone(target);
      tri::UpdateSelection<MyMesh>::FaceFromBorderFlag(target);
    }
    Rvcg::KDtree< MyMesh, PcMesh >::getBary(target, bary);
    List indices = Rvcg::KDtree< PcMesh, MyMesh >::KDtreeIO(bary, query, k,nofP, mDepth,threads);
    IntegerMatrix ktree = indices["index"];
    NumericMatrix iomat(3,query.vn), normals(3,query.vn);
    NumericMatrix barycoord(3,query.vn);
    IntegerVector border(query.vn), faceptr(query.vn);
    std::fill(border.begin(), border.end(),0);
    //MyMesh::VertexIterator vi = query.vert.begin();
    NumericVector distout(query.vn);
    NumericVector angle(query.vn);
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
#pragma omp parallel for schedule(static)
    for (int i = 0; i < query.vn; i++) {
      MyMesh::VertexIterator vi = query.vert.begin()+i;
      Point3f clost;
      MyMesh::CoordType tt, tmpnorm;
      MyFace::ScalarType dist = 1e12;
      MyFace::ScalarType distance_old = 1e12;
      Point3f currp = (*vi).P();
      clost = currp;
      float ang;
      Point3f tmp = (*vi).P();
      for (int j=0; j < k; j++) {
	if (ktree(i,j) != -1) {
	  dist = 1e12;
	  int fptr = ktree(i,j);
	  PointDistanceBase(target.face[fptr],currp, dist, tmp);
	  
	  // check normal deviation and if it deviates from threshold, hit point is discarded
	  // if no point matches criterion, the closest point on the face with the closest barycenter is selected
	  // and distance set to 1e5
	  if (  0 < angdev) { 
	    
	    MyMesh::CoordType refnorm = (*vi).N();
	    tmpnorm = clost*0;
	    if (!wnorm) {
	      tmpnorm = target.face[fptr].N();
	    } else {
	      for (int j1=0; j1 <3;j1++) {
		Point3f vdist = target.face[fptr].V(j1)->P() - tmp;
		float weight = sqrt(vdist.dot(vdist));
		if (weight > 0)
		  weight = 1/weight;
		else 
		  weight = 1e12;
		tmpnorm += target.face[fptr].V(j1)->N()*weight;
	      }
	    }
	    ang = Angle(tmpnorm,refnorm);
	    
	    if (ang > angdev)
	      dist = 1e5;
	  }
	  if (dist < distance_old) {
	    distance_old = dist;
	    clost = tmp;
	    faceptr[i] = fptr;
	    if (angdev > 0)
	      angle[i] = ang;
	  }
	}
      }
      distout[i] = distance_old;
      tt = clost*0;
	    
      // get normals at hit point
      if (facenormals) {
	tt = target.face[faceptr[i]].N();
      } else {
	for (int j=0; j <3;j++) {
	  //if (&(target.face[faceptr[i]].V(j)->N())) {
	  Point3f vdist = target.face[faceptr[i]].V(j)->P() - clost;
	  float weight = sqrt(vdist.dot(vdist));
	  if (weight > 0)
	    weight = 1/weight;
	  else 
	    weight = 1e12;
	    
	  tt +=(target.face[faceptr[i]].V(j)->N()*weight);
	}
	//}
      }
      float vl = sqrt(tt.dot(tt));
      if (vl > 0) {//check for zero length normals
	tt=tt/vl;
      }   
      // calculate sign for distances
      if (sign && (distout[i] < 1e5)) {
	Point3f dif = clost - currp;
	//float sign = dif.dot(tt);	
	if (dif.dot(tt) < 0)
	  distout[i] = -distout[i] ;
      }
      // write back
      iomat(0,i) = clost[0];
      iomat(1,i) = clost[1];
      iomat(2,i) = clost[2];
      normals(0,i) = tt[0];
      normals(1,i) = tt[1];    
      normals(2,i) = tt[2];
      // get barycentric coordinates
      if (barycentric) {
      
	MyMesh::CoordType baryco = currp*0;
	InterpolationParameters<MyFace,ScalarType>(target.face[faceptr[i]],target.face[faceptr[i]].N(),clost,baryco);
	barycoord(0,i) = baryco[0];
	barycoord(1,i) = baryco[1];
	barycoord(2,i) = baryco[2];
      }
      if (borderchk) {
	if (target.face[faceptr[i]].IsS())
	  border[i] = 1;
      }
     
    }
    return List::create(Named("iomat")=iomat,
			Named("distance")= distout,
			Named("faceptr")= faceptr,
			Named("barycoord") = barycoord,
			Named("border") = border, 
			Named("normals") = normals,
			Named("angle") = angle
			);
      
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP Rbarycenter(SEXP vb_, SEXP it_) {
  try {
    MyMesh m;
    int checkit = Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    MyMesh out;
    Rvcg::KDtree< MyMesh, MyMesh >::getBary(m, out);
    Rcpp::NumericMatrix barycoord(3,out.vn);
    for (int i = 0; i < out.vn; i++) {
      MyMesh::CoordType tmp;
      tmp = out.vert[i].cP();
      barycoord(0,i) = tmp[0];
      barycoord(1,i) = tmp[1];
      barycoord(2,i) = tmp[2];
    }
    return wrap(barycoord);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
