
#ifndef TDHAUSDORF
#define TDHAUSDORF

#include <time.h>


#include <vcg/math/histogram.h>
#include <vcg/complex/complex.h>
#include <vcg/simplex/face/component_ep.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/complex/algorithms/update/component_ep.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include "metroSampling.h"

using namespace std;
// #include "/media/celli/Data/Rvcg/src/vcglib/apps/metro/sampling.h"
// using namespace std;

// project definitions.
// error messages

#define MSG_ERR_MESH_LOAD               "error loading the input meshes.\n"
#define MSG_ERR_INVALID_OPTION          "unable to parse option '%s'\n"
#define MSG_ERR_FILE_OPEN               "unable to open the output file.'n"
#define MSG_ERR_UNKNOWN_FORMAT          "unknown file format '%s'.\n"

// global constants
#define NO_SAMPLES_PER_FACE             10
#define N_SAMPLES_EDGE_TO_FACE_RATIO    0.1
#define BBOX_FACTOR                     0.1
#define INFLATE_PERCENTAGE			    0.02
#define MIN_SIZE					    125		/* 125 = 5^3 */
#define N_HIST_BINS                     256
#define PRINT_EVERY_N_ELEMENTS          1000


class CFaceMetro;
class CVertexMetro;
struct UsedTypes:public vcg::UsedTypes< vcg::Use<CFaceMetro>::AsFaceType, vcg::Use<CVertexMetro>::AsVertexType>{};
class CVertexMetro   : public vcg::Vertex<UsedTypes,vcg::vertex::Coord3d,vcg::vertex::Qualityf,vcg::vertex::Normal3d,vcg::vertex::Color4b,vcg::vertex::BitFlags> {};
class CFaceMetro     : public vcg::Face< UsedTypes,vcg::face::VertexRef, vcg::face::Normal3d, vcg::face::EdgePlane,vcg::face::Color4b,vcg::face::Mark,vcg::face::BitFlags> {};
class CMeshMetro     : public vcg::tri::TriMesh< std::vector<CVertexMetro>, std::vector<CFaceMetro> > {};

#endif
