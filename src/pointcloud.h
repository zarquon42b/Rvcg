#include<vcg/complex/complex.h>
using namespace vcg;
//using namespace std;
// The class prototypes.
//we are setting up a mesh without edges or borders, just plain pointclouds

class PcVertex; 
struct PcUsedTypes: public UsedTypes<Use<PcVertex>::AsVertexType>{};

class PcEdge : public Edge<PcUsedTypes>{};
class PcVertex  : public Vertex< PcUsedTypes, 
		  //vertex::InfoOcf,
                                 vertex::Coord3f, 
                                 vertex::BitFlags, 
                                 vertex::Normal3f, 
                                 vertex::Mark,
                                 vertex::Color4b, 
                                 vertex::Qualityf >{};


// ocf class
//class PcMesh : public vcg::tri::TriMesh< vcg::vertex::vector_ocf<PcVertex>, vcg::face::vector_ocf<PcFace > >{};
// default class
class PcMesh : public vcg::tri::TriMesh< std::vector<PcVertex> >{};
typedef  PcMesh::ScalarType PCScalarType;
typedef  PcMesh::VertexIterator PCVertexIterator;
typedef  PcMesh::VertexPointer PCVertexPointer;
typedef  PcMesh::CoordType PCCoordType;
typedef  PcMesh::ScalarType PCScalarType;
typedef  PcMesh::VertContainer PCVertContainer;

/*typedef PcMesh::ConstVertexIterator ConstVertexIterator;
  typedef PcMesh::ConstFaceIterator   ConstFaceIterator;*/
