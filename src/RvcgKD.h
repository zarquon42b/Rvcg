#include <RcppArmadillo.h>
#include <vector>
#include <vcg/complex/complex.h>
#include <vcg/space/index/kdtree/kdtree.h>
#include <Rconfig.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace tri;
using namespace Rcpp;


namespace Rvcg
{
  template <class MeshTarget, class MeshQuery>
  class KDtree
  {
  public:
    typedef typename MeshTarget::CoordType      CoordTarget;
    typedef typename MeshTarget::ScalarType     ScalarTarget;
    typedef typename MeshTarget::VertexType     VertexTarget;
    typedef typename MeshTarget::VertexIterator VertexIteratorTarget;
    typedef typename MeshTarget::FaceType       FaceTarget;
    typedef typename MeshTarget::FacePointer    FacePointerTarget;
    typedef typename MeshTarget::FaceIterator   FaceIteratorTarget;
    typedef typename MeshTarget::FaceContainer  FaceContainerTarget;
    typedef typename MeshTarget::VertContainer  VertContainerTarget;


    typedef typename MeshQuery::CoordType      CoordQuery;
    typedef typename MeshQuery::ScalarType     ScalarQuery;
    typedef typename MeshQuery::VertexPointer  VertexQuery;
    typedef typename MeshQuery::VertexIterator VertexIteratorQuery;
    typedef typename MeshQuery::FaceType       FaceQuery;
    typedef typename MeshQuery::FacePointer    FacePointerQuery;
    typedef typename MeshQuery::FaceIterator   FaceIteratorQuery;
    typedef typename MeshQuery::FaceContainer  FaceContainerQuery;

    static KdTree<float> KDTreeCreate(MeshTarget &target, unsigned int nofPointsPerCell, unsigned int maxDepth) {
      try {
	VertexConstDataWrapper<MeshTarget> ww(target);
	KdTree<float> tree(ww, nofPointsPerCell, maxDepth);
	return tree;
      } catch (std::exception& e) {
	::Rf_error( e.what());
      } catch (...) {
	::Rf_error("unknown exception");
      }
    }
    static  List KDtreeIO(MeshTarget &target, MeshQuery &query, int k, unsigned int nofPointsPerCell = 16, unsigned int maxDepth = 64, int threads = 1) {
      try {
	typedef pair<float,int> mypair;
	IntegerMatrix result(query.vn,k);
	NumericMatrix distance(query.vn,k);
	std::fill(result.begin(), result.end(),-1);
	KdTree<float> tree = KDTreeCreate(target, nofPointsPerCell, maxDepth);
	//tree.setMaxNofNeighbors(k);
	KdTree<float>::PriorityQueue queue;
#pragma omp parallel for firstprivate(queue, tree) schedule(static) num_threads(threads)
	for (int i = 0; i < query.vn; i++) {
	  //tree.doQueryK(query.vert[i].cP());
	  tree.doQueryK(query.vert[i].cP(), k, queue);
	  //int neighbours = tree.getNofFoundNeighbors();
	  int neighbours = queue.getNofElements();
	  vector<mypair> sortit;
	  for (int j=0; j < neighbours; j++) {      
	    int neightId = queue.getIndex(j);
	    float dist = Distance(query.vert[i].cP(),target.vert[neightId].cP());
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
    
    static void getBary(MeshTarget &inmesh, MeshQuery &outmesh) {
      vcg::tri::Allocator<MeshQuery>::AddVertices(outmesh,inmesh.fn);
      VertexIteratorQuery vi = outmesh.vert.begin();
      for (int i = 0; i < inmesh.fn; i++) {
	(*vi).P() =  Barycenter<FaceTarget>(inmesh.face[i]);
	++vi;
      }
    }

    static List clostKD(MeshTarget &target, MeshQuery &query, arma::imat &closest_indices, int k, double angdev, bool facenormals, bool sign, bool weightnorm , bool borderchk, bool barycentric,int threads) {
      try {
	//setup output datastructures
	arma::mat ioclost(4,query.vn), normals(4,query.vn);
	ioclost.fill(1);normals.fill(1);
	arma::mat barycoord;
	arma::ivec border, faceptr(query.vn);
	arma::vec distances(query.vn);
	float mynan = std::nan("1");
	if (barycentric)
	  barycoord.resize(3,query.vn);
	if (borderchk) {
	  border.resize(query.vn);
	  border.fill(0);
	}
	arma::vec angle;
	if (angdev > 0)
	  angle.resize(query.vn);
#pragma omp parallel for schedule(static) num_threads(threads)
	for (int i = 0; i < query.vn; i++) {
	  MyMesh::VertexIterator vi = query.vert.begin()+i;
	  Point3f clost;
	  MyMesh::CoordType vertexnormal, tmp_vertexnormal;
	  MyFace::ScalarType dist = 1e12;
	  MyFace::ScalarType distance_old = 1e12;
	  Point3f currp = (*vi).P();
	  clost = currp;
	  float ang;
	  Point3f tmp = (*vi).P();
	  for (int j=0; j < k; j++) {
	    if (closest_indices(i,j) != -1) {
	      dist = 1e12;
	      int fptr = closest_indices(i,j);
	      PointDistanceBase(target.face[fptr],currp, dist, tmp);
	  
	      // check normal deviation and if it deviates from threshold, hit point is discarded
	      // if no point matches criterion, the closest point on the face with the closest barycenter is selected
	      // and distance set to 1e5
	      if (  0 < angdev) { 
	    
		MyMesh::CoordType refnorm = (*vi).N();
		tmp_vertexnormal = clost*0;
		if (!weightnorm) {
		  tmp_vertexnormal = target.face[fptr].N();
		} else {
		  for (int j1=0; j1 <3;j1++) {
		    Point3f vdist = target.face[fptr].V(j1)->P() - tmp;
		    float weight = sqrt(vdist.dot(vdist));
		    if (weight > 0)
		      weight = 1/weight;
		    else 
		      weight = 1e12;
		    tmp_vertexnormal += target.face[fptr].V(j1)->N()*weight;
		  }
		}
		ang = Angle(tmp_vertexnormal,refnorm);
	    
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
	  distances[i] = distance_old;
	  vertexnormal = clost*0;
	    
	  // get normals at hit point
	  if (facenormals) {//use face normals
	    vertexnormal = target.face[faceptr[i]].N();
	  } else { //use weighted vertex normals
	    for (int j=0; j <3;j++) {
	      Point3f vdist = target.face[faceptr[i]].V(j)->P() - clost;
	      float weight = sqrt(vdist.dot(vdist));
	      if (weight > 0)
		weight = 1/weight;
	      else 
		weight = 1e12;
	      vertexnormal += target.face[faceptr[i]].V(j)->N()*weight;
	    }
	  }
      
	  float vl = sqrt(vertexnormal.dot(vertexnormal));
	  if (vl > 0) {//check for zero length normals
	    vertexnormal=vertexnormal/vl;
	  }   
	  // calculate sign for distances
	  if (sign && (distances[i] < 1e12)) {
	    Point3f dif = clost - currp;
	    //float sign = dif.dot(vertexnormal);	
	    if (dif.dot(vertexnormal) < 0)
	      distances[i] = -distances[i];
	  } else if (distances[i] >= 1e12) {
	    distances[i] = mynan;
	  }
	  // write back
	  for (int j=0; j < 3;j++) {
	    ioclost(j,i) = clost[j];
	    normals(j,i) = vertexnormal[j];
	  }
	  // get barycentric coordinates
	  if (barycentric) {
      	    MyMesh::CoordType baryco = currp*0;
	    InterpolationParameters<MyFace,ScalarType>(target.face[faceptr[i]],target.face[faceptr[i]].N(),clost,baryco);
	    for (int j=0; j < 3;j++)
	      barycoord(j,i) = baryco[j];
	  }
	  if (borderchk) {
	    if (target.face[faceptr[i]].IsS())
	      border[i] = 1;
	  }
     	}
	faceptr = faceptr+1;
	
	List out = List::create(Named("vb")=ioclost,
				Named("it")=wrap(1),
				Named("normals") = normals,
				Named("quality")= NumericVector(distances.begin(),distances.end()),
				Named("faceptr")= NumericVector(faceptr.begin(),faceptr.end())
				);
			    
	if (barycentric)
	  out["barycoords"] = barycoord;
	if (borderchk)
	  out["border"] = NumericVector(border.begin(),border.end());
	if (angdev > 0)
	  out["angle"] = NumericVector(angle.begin(),angle.end());
	out.attr("class") = "mesh3d";
	return out;
			   
      } catch (std::exception& e) {
	::Rf_error( e.what());
      } catch (...) {
	::Rf_error("unknown exception");
      }
    }
  };
}
