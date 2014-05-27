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
class CurveMyFace;
class CurveMyEdge;
class CurveMyVertex;
struct CurveMyUsedTypes: public UsedTypes<Use<CurveMyVertex>::AsVertexType,
  Use<CurveMyEdge>::AsEdgeType,
  Use<CurveMyFace>::AsFaceType
  >{};

class CurveMyEdge : public Edge<CurveMyUsedTypes>{};
class CurveMyVertex  : public Vertex< CurveMyUsedTypes, 
				      //vertex::InfoOcf,
				      vertex::Coord3f, 
				      vertex::BitFlags, 
				      vertex::Normal3f, 
				      vertex::Mark,
				      //vertex::Color4b, 
				      vertex::Qualityf,
				      vertex::VFAdj,
				      vertex::Curvaturef,
				      vertex::CurvatureDirf
				 >{};
class CurveMyFace: public Face  <CurveMyUsedTypes, 
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
//class CurveMyMesh : public vcg::tri::TriMesh< vcg::vertex::vector_ocf<CurveMyVertex>, vcg::face::vector_ocf<CurveMyFace > >{};
// default class
class CurveMyMesh : public vcg::tri::TriMesh< std::vector<CurveMyVertex>, std::vector<CurveMyFace > >{};
typedef  CurveMyMesh::ScalarType ScalarType;
typedef  CurveMyMesh::VertexIterator VertexIterator;
typedef  CurveMyMesh::VertexPointer VertexPointer;
typedef  CurveMyMesh::FaceIterator FaceIterator;
typedef  CurveMyMesh::FacePointer FacePointer;
typedef  CurveMyMesh::EdgePointer   EdgePointer;
typedef  CurveMyMesh::EdgeIterator   EdgeIterator;
typedef  CurveMyMesh::CoordType CoordType;
typedef  CurveMyMesh::ScalarType ScalarType;
typedef  CurveMyMesh::VertContainer VertContainer;
typedef  CurveMyMesh::FaceContainer FaceContainer;
/*typedef CurveMyMesh::ConstVertexIterator ConstVertexIterator;
  typedef MyMesh::ConstFaceIterator   ConstFaceIterator;*/
