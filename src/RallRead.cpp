#include "typedefImport.h"
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/container/simple_temporary_data.h>
#include <vcg/complex/algorithms/update/texture.h>
#include <vcg/complex/algorithms/attribute_seam.h>
#include <string.h>
#include <Rcpp.h>  

using namespace Rcpp;

/// The following to helper functions are copied from filter_texture plugin of meshlab
/////// FUNCTIONS NEEDED BY "UV WEDGE TO VERTEX" FILTER
inline void ExtractVertex(const MyMeshImport & srcMesh, const MyMeshImport::FaceType & f, int whichWedge, const MyMeshImport & dstMesh, MyMeshImport::VertexType & v)
{
    (void)srcMesh;
    (void)dstMesh;
    // This is done to preserve every single perVertex property
    // perVextex Texture Coordinate is instead obtained from perWedge one.
    v.ImportData(*f.cV(whichWedge));
    v.T() = f.cWT(whichWedge);
}

inline bool CompareVertex(const MyMeshImport & m, const MyMeshImport::VertexType & vA, const MyMeshImport::VertexType & vB)
{
    (void)m;
    return (vA.cT() == vB.cT());
}
///////

RcppExport SEXP RallRead(SEXP filename_, SEXP updateNormals_, SEXP colorread_, SEXP clean_) 
{
  std::string str = Rcpp::as<std::string>(filename_);
  const char *filename = str.c_str();
  bool updateNormals = as<bool>(updateNormals_);
  bool colorread = as<bool>(colorread_);
  bool clean = as<bool>(clean_);
  MyMeshImport m;
  int err2 = tri::io::Importer<MyMeshImport>::Open(m,filename);
  if (err2) {
    return wrap(1);
  } else { 
    if (m.fn == 0)
      updateNormals = false;
    SimpleTempData<MyMeshImport::VertContainer,int> indices(m.vert);
   
    if (clean) {
      int dup = tri::Clean<MyMeshImport>::RemoveDuplicateVertex(m);
      int dupface = tri::Clean<MyMeshImport>::RemoveDuplicateFace(m);
      int unref =  tri::Clean<MyMeshImport>::RemoveUnreferencedVertex(m);
      vcg::tri::Allocator< MyMeshImport >::CompactVertexVector(m);
      vcg::tri::Allocator< MyMeshImport >::CompactFaceVector(m);
      if (dup > 0 || unref > 0 || dupface > 0)
	Rprintf("Removed %i duplicate %i unreferenced vertices and %i duplicate faces\n",dup,unref,dupface);
    }  
    // do texture processing
    bool tex = false;
    std::vector<float> texvec;
    std::vector<string> texfile;
    if (m.textures.size() > 0) {
      tex = true;
      tri::AttributeSeam::SplitVertex(m, ExtractVertex, CompareVertex);
      texfile = m.textures;
      texvec.resize(2*m.vn);
      vcg::tri::Allocator< MyMeshImport >::CompactVertexVector(m);
      vcg::tri::Allocator< MyMeshImport >::CompactFaceVector(m);
    }
   
    // setup output structures
    NumericVector vb(3*m.vn);    
    std::vector<int> colvec;
    if (colorread)
      colvec.resize(3*m.vn);
    
    IntegerVector it(3*m.fn);
    std::vector<double> normals;
    if (updateNormals)
      normals.resize(3*m.vn);
    
    if (updateNormals) { // update Normals
      tri::UpdateNormal<MyMeshImport>::PerVertexNormalized(m);
    }
    // write back
    VertexIterator vi=m.vert.begin();
    for (int i=0;  i < m.vn; i++) {
      vb(i*3) = (*vi).P()[0];
      vb(i*3+1) = (*vi).P()[1];
      vb(i*3+2) = (*vi).P()[2];
      indices[vi] = i;
      if (updateNormals) {
	normals[i*3] = (*vi).N()[0];
	normals[i*3+1] = (*vi).N()[1];
	normals[i*3+2] = (*vi).N()[2];
      }
      if (colorread) {
	colvec[i*3] = (*vi).C()[0];
	colvec[i*3+1] = (*vi).C()[1];
	colvec[i*3+2] = (*vi).C()[2];
      }
      if (tex) {
	texvec[i*2] = (*vi).T().U();
	texvec[i*2+1] = (*vi).T().V();
      }
      ++vi;
    }
    FacePointer fp;
    int vv[3];
    FaceIterator fi=m.face.begin();
    if (m.fn > 0) {
      for (int i=0; i < m.fn;i++) {
	fp=&(*fi);
	if( ! fp->IsD() ) {
	  vv[0]=indices[fp->cV(0)];
	  vv[1]=indices[fp->cV(1)];
	  vv[2]=indices[fp->cV(2)];
	  it(i*3)=vv[0];
	  it(i*3+1)=vv[1];
	  it(i*3+2)=vv[2];
	  ++fi;
	}
      }
    }
    //return wrap(vb);
    return List::create(Named("vb") = vb, 
			Named("it") = it,
			Named("normals") = normals,
			Named("colors") = colvec,
			Named("texcoord") = texvec,
			Named("texfile") = texfile
			);
  }
 }
