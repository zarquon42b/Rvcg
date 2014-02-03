#include <Rcpp.h>
#include <vector>
#include <vcg/complex/complex.h>
#include <vcg/space/index/kdtree/kdtree.h>

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

    static  List KDtreeIO(MeshTarget &target, MeshQuery &query, int k, unsigned int nofPointsPerCell = 16, unsigned int maxDepth = 64) {
      typedef pair<float,int> mypair;
      IntegerMatrix result(query.vn,k);
      NumericMatrix distance(query.vn,k);
      std::fill(result.begin(), result.end(),-1);
      VertexConstDataWrapper<MeshTarget> ww(target);
      KdTree<float> tree(ww, nofPointsPerCell, maxDepth);
      tree.setMaxNofNeighbors(k);
      for (int i = 0; i < query.vn; i++) {
	tree.doQueryK(query.vert[i].cP());
	int neighbours = tree.getNofFoundNeighbors();
	vector<mypair> sortit;
	for (int j=0; j < neighbours; j++) {      
	  int neightId = tree.getNeighborId(j);
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
    }
    
    static void getBary(MeshTarget &inmesh, MeshQuery &outmesh) {
      vcg::tri::Allocator<MeshQuery>::AddVertices(outmesh,inmesh.fn);
      VertexIteratorQuery vi = outmesh.vert.begin();
      for (int i = 0; i < inmesh.fn; i++) {
	(*vi).P() =  Barycenter<FaceTarget>(inmesh.face[i]);
	++vi;
      }
    }
    
    
  };
}
