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
#include <wrap/io_trimesh/export_vrml.h>

#include <vcg/complex/algorithms/update/color.h>
//#include <vcg/complex/algorithms/update/curvature.h>

#ifndef RcppExport
#define RcppExport extern "C"
#endif

using namespace vcg;
//using namespace std;
// The class prototypes.
class MyFaceImport;
class MyEdgeImport;
class MyVertexImport;
struct MyUsedTypes: public UsedTypes<Use<MyVertexImport>::AsVertexType,
                                     Use<MyEdgeImport>::AsEdgeType,
                                     Use<MyFaceImport>::AsFaceType
                                     >{};

class MyEdgeImport : public Edge<MyUsedTypes>{};
class MyVertexImport  : public Vertex< MyUsedTypes, 
				       vertex::InfoOcf,
				       vertex::Coord3f, 
				       vertex::BitFlags, 
				       vertex::Normal3fOcf, 
				       vertex::Mark,
				       vertex::Color4bOcf, 
				       vertex::QualityfOcf,
				       /* vertex::VFAdjOcf,
					  vertex::CurvaturefOcf,
					  vertex::CurvatureDirfOcf,*/
				       vertex::VFAdjOcf,
				       vertex::TexCoordfOcf
				       >{};
class MyFaceImport: public Face  <MyUsedTypes, 
				  face::InfoOcf,
				  //	  face::VFAdj,
				  face::VertexRef,
				  face::BitFlags,
				  face::Mark,
				  //face::FFAdjOcf, 
				  face::WedgeTexCoordfOcf,
				  face::Color4bOcf, 
				  // face::FFAdj, 
				  face::Normal3fOcf,
				  face::QualityfOcf
				  >{};


// ocf class
class MyMeshImport : public vcg::tri::TriMesh< vcg::vertex::vector_ocf<MyVertexImport>, vcg::face::vector_ocf<MyFaceImport > >{};
// default class
//class MyMeshImport : public vcg::tri::TriMesh< std::vector<MyVertexImport>, std::vector<MyFaceImport > >{};
typedef  MyMeshImport::ScalarType ScalarType;
typedef  MyMeshImport::VertexIterator VertexIterator;
typedef  MyMeshImport::VertexPointer VertexPointer;
typedef  MyMeshImport::FaceIterator FaceIterator;
typedef  MyMeshImport::FacePointer FacePointer;
typedef  MyMeshImport::EdgePointer   EdgePointer;
typedef  MyMeshImport::EdgeIterator   EdgeIterator;
typedef  MyMeshImport::CoordType CoordType;
typedef  MyMeshImport::ScalarType ScalarType;
typedef  MyMeshImport::VertContainer VertContainer;
typedef  MyMeshImport::FaceContainer FaceContainer;
/*typedef MyMesh::ConstVertexIterator ConstVertexIterator;
  typedef MyMesh::ConstFaceIterator   ConstFaceIterator;*/
