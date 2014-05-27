#include <string.h>
#include <vector>
#include <stdio.h>
#include <cstddef>
//#include <iostream>

// VCG headers for triangular mesh processing

#include <vcg/complex/complex.h>
// include update algorithms
#include <vcg/complex/algorithms/update/topology.h>
//#include <vcg/complex/algorithms/update/edges.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/update/normal.h>

#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/intersection.h>
//include headers for search grids
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/smooth.h>
//#include <vcg/complex/allocate.h>
#include <wrap/callback.h>
#include <vcg/complex/append.h>
#include <vcg/container/simple_temporary_data.h>

// VCG File Format Importer/Exporter
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/update/color.h>
//#include <vcg/complex/algorithms/update/curvature.h>
#ifndef RcppExport
#define RcppExport extern "C"
#endif

using namespace vcg;
//using namespace std;
// The class prototypes.
class TopoMyFace;
class TopoMyEdge;
class TopoMyVertex;
struct TopoMyUsedTypes: public UsedTypes<Use<TopoMyVertex>::AsVertexType,
  Use<TopoMyEdge>::AsEdgeType,
  Use<TopoMyFace>::AsFaceType
  >{};

class TopoMyEdge : public Edge<TopoMyUsedTypes>{};
class TopoMyVertex  : public Vertex< TopoMyUsedTypes, 
				 //vertex::InfoOcf,
				 vertex::Coord3f, 
				 vertex::BitFlags, 
				 vertex::Normal3f, 
				 vertex::Mark,
				 //vertex::Color4b, 
				 //vertex::Qualityf,
				 vertex::VFAdj
				 //vertex::Curvaturef,
				 //vertex::CurvatureDirf
				 >{};
class TopoMyFace: public Face  <TopoMyUsedTypes, 
			    //face::InfoOcf,
			    face::VertexRef,
			    face::BitFlags,
			    face::Mark,
			    face::FFAdj, 
			    face::VFAdj, 
			    face::Normal3f
			    /*face::FFAdj, 
			      face::VFAdj,*/
			    > {};


// ocf class
//class TopoMyMesh : public vcg::tri::TriMesh< vcg::vertex::vector_ocf<TopoMyVertex>, vcg::face::vector_ocf<TopoMyFace > >{};
// default class
class TopoMyMesh : public vcg::tri::TriMesh< std::vector<TopoMyVertex>, std::vector<TopoMyFace > >{};
typedef  TopoMyMesh::ScalarType ScalarType;
typedef  TopoMyMesh::VertexIterator VertexIterator;
typedef  TopoMyMesh::VertexPointer VertexPointer;
typedef  TopoMyMesh::FaceIterator FaceIterator;
typedef  TopoMyMesh::FacePointer FacePointer;
typedef  TopoMyMesh::EdgePointer   EdgePointer;
typedef  TopoMyMesh::EdgeIterator   EdgeIterator;
typedef  TopoMyMesh::CoordType CoordType;
typedef  TopoMyMesh::ScalarType ScalarType;
typedef  TopoMyMesh::VertContainer VertContainer;
typedef  TopoMyMesh::FaceContainer FaceContainer;
/*typedef TopoMyMesh::ConstVertexIterator ConstVertexIterator;
  typedef MyMesh::ConstFaceIterator   ConstFaceIterator;*/
