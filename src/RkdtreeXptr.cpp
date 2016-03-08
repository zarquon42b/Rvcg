#include "typedef.h"
#include "pointcloud.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
#include <RvcgKD.h>

RcppExport SEXP createKDtree(SEXP target_, SEXP nofPointsPerCell_, SEXP maxDepth_) {
  Rcpp::XPtr< MyMesh > target(new MyMesh,true);
  Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(*target,target_);
  unsigned int nofPointsPerCell = as<unsigned int >(nofPointsPerCell_);
  unsigned int maxDepth = as<unsigned int >(maxDepth_);
  VertexConstDataWrapper<MyMesh> ww(*target);
  KdTree<float> tree(ww, nofPointsPerCell, maxDepth);

  Rcpp::XPtr< KdTree<float> > Rtree(new KdTree<float>(ww, nofPointsPerCell, maxDepth),true);
  
  
  return List::create(Named("kdtree") =Rtree,
		      Named("target") = target);

}

List searchKDtree(XPtr< KdTree<float> > kdtree, XPtr< MyMesh > target, MyMesh &query, int k, int threads){
  try {
    KdTree<float>::PriorityQueue queue;
    typedef pair<float,int> mypair;
    IntegerMatrix result(query.vn,k);
    NumericMatrix distance(query.vn,k);
    std::fill(result.begin(), result.end(),-1);
#pragma omp parallel for firstprivate(queue, kdtree) schedule(static) num_threads(threads)
    for (int i = 0; i < query.vn; i++) {
      //tree.doQueryK(query.vert[i].cP());
      (*kdtree).doQueryK(query.vert[i].cP(), k, queue);
      //int neighbours = tree.getNofFoundNeighbors();
      int neighbours = queue.getNofElements();
      vector<mypair> sortit;
      for (int j=0; j < neighbours; j++) {      
	int neightId = queue.getIndex(j);
	float dist = Distance(query.vert[i].cP(),(*target).vert[neightId].cP());
	sortit.push_back(mypair(dist, neightId));
      }

      sort(sortit.begin(),sortit.end());
      for (int j = 0; j < neighbours; j++){
	result(i,j) = sortit[j].second;
	distance(i,j) = sortit[j].first;
      }
    }
    return List::create(Named("index") = result,
			Named("distance") = distance)
      ;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP RsearchKDtree(SEXP kdtree_,SEXP target_, SEXP query_, SEXP k_, SEXP threads_) {
  try {
    XPtr< KdTree<float> > kdtree(kdtree_);
    XPtr< MyMesh > target(target_);
    MyMesh query;
    Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(query,query_);
    int k = as<int>(k_);
    int threads = as<int>(threads_);
    List out = searchKDtree(kdtree, target, query, k, threads);
    return out;
    // start searching
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP searchKDtreeForClosestPoints(SEXP kdtree_,SEXP bary_, SEXP targetmesh_, SEXP query_, SEXP k_,SEXP sign_,SEXP borderchk_,SEXP barycentric_, SEXP angdev_=wrap(0), SEXP wnorm_=wrap(true), SEXP facenormals_=wrap(false), SEXP threads_=wrap(1)) {
  try {
    XPtr< KdTree<float> > kdtree(kdtree_);
    XPtr< MyMesh > bary(bary_);
    XPtr< MyMesh > target(targetmesh_);
    MyMesh query;
    Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(query,query_);
    int k = as<int>(k_);
    int threads = as<int>(threads_);
    bool barycentric = as<bool>(barycentric_);
    bool borderchk = as<bool>(borderchk_);
    bool wnorm = as<bool>(wnorm_);
    double angdev = as<double>(angdev_);
    bool sign = as<bool>(sign_);
    bool facenormals = as<bool>(facenormals_);
    (*target).face.EnableNormal();
    if (angdev > 0) {
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(query);
    }
    tri::UpdateNormal<MyMesh>::PerFaceNormalized(*target);
    tri::UpdateNormal<MyMesh>::PerVertexNormalized(*target);
    if (borderchk) { //update Border flags
      tri::UpdateFlags<MyMesh>::FaceBorderFromNone(*target);
      tri::UpdateSelection<MyMesh>::FaceFromBorderFlag(*target);
    }
    List indices = searchKDtree(kdtree,bary,query,k,threads);
    IntegerMatrix ktree = indices["index"];
    NumericMatrix iomat(3,query.vn), normals(3,query.vn);
    NumericMatrix barycoord(3,query.vn);
    IntegerVector border(query.vn), faceptr(query.vn);
    std::fill(border.begin(), border.end(),0);
    //MyMesh::VertexIterator vi = query.vert.begin();
    NumericVector distout(query.vn);
    NumericVector angle(query.vn);
#pragma omp parallel for schedule(static) num_threads(threads)
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
	  PointDistanceBase((*target).face[fptr],currp, dist, tmp);
	  
	  // check normal deviation and if it deviates from threshold, hit point is discarded
	  // if no point matches criterion, the closest point on the face with the closest barycenter is selected
	  // and distance set to 1e5
	  if (  0 < angdev) { 
	    
	    MyMesh::CoordType refnorm = (*vi).N();
	    tmpnorm = clost*0;
	    if (!wnorm) {
	      tmpnorm = (*target).face[fptr].N();
	    } else {
	      for (int j1=0; j1 <3;j1++) {
		Point3f vdist = (*target).face[fptr].V(j1)->P() - tmp;
		float weight = sqrt(vdist.dot(vdist));
		if (weight > 0)
		  weight = 1/weight;
		else 
		  weight = 1e12;
		tmpnorm += (*target).face[fptr].V(j1)->N()*weight;
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
       tt = (*target).face[faceptr[i]].N();
     } else {
       for (int j=0; j <3;j++) {
	 //if (&((*target).face[faceptr[i]].V(j)->N())) {
	 Point3f vdist = (*target).face[faceptr[i]].V(j)->P() - clost;
	 float weight = sqrt(vdist.dot(vdist));
	 if (weight > 0)
	   weight = 1/weight;
	 else 
	   weight = 1e12;
	 
	 tt +=((*target).face[faceptr[i]].V(j)->N()*weight);
       }
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
	InterpolationParameters<MyFace,ScalarType>((*target).face[faceptr[i]],(*target).face[faceptr[i]].N(),clost,baryco);
	barycoord(0,i) = baryco[0];
	barycoord(1,i) = baryco[1];
	barycoord(2,i) = baryco[2];
      }
      if (borderchk) {
	if ((*target).face[faceptr[i]].IsS())
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
  


