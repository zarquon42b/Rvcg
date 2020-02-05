#include "typedef.h"

#include <RcppArmadillo.h>  
using namespace Rcpp;

/// The following to helper functions are copied from filter_texture plugin of meshlab
/////// FUNCTIONS NEEDED BY "UV WEDGE TO VERTEX" FILTER
void ExtractVertex(const MyMesh & srcMesh, const MyMesh::FaceType & f, int whichWedge, const MyMesh & dstMesh, MyMesh::VertexType & v)
{
  (void)srcMesh;
  (void)dstMesh;
  // This is done to preserve every single perVertex property
  // perVextex Texture Coordinate is instead obtained from perWedge one.
  v.ImportData(*f.cV(whichWedge));
  v.T() = f.cWT(whichWedge);
}

bool CompareVertex(const MyMesh & m, const MyMesh::VertexType & vA, const MyMesh::VertexType & vB)
{
  (void)m;
  return (vA.cT() == vB.cT());
}
///////

RcppExport SEXP RallRead(SEXP filename_, SEXP updateNormals_, SEXP colorread_, SEXP clean_,SEXP silent_, SEXP type_) 
{
  try {
    std::string str = Rcpp::as<std::string>(filename_);
    const char *filename = str.c_str();
    bool updateNormals = as<bool>(updateNormals_);
    bool colorread = as<bool>(colorread_);
    bool clean = as<bool>(clean_);
    bool silent = as<bool>(silent_);
    int type = as<int>(type_);
    MyMesh m; 
    bool hasNormal = false, WedgeTex=false, VertTex = false, hasQuality=false,hasFaceQuality=false;
    int mask0 = 0; //initializie import mask
    tri::io::Importer<MyMesh>::LoadMask(filename, mask0);
    // start allocating space for availables stuff
    if( (mask0 & tri::io::Mask::IOM_VERTCOLOR) && colorread)
      m.vert.EnableColor();
    if( (mask0 & tri::io::Mask::IOM_WEDGTEXCOORD && colorread)) {
      m.face.EnableWedgeTexCoord();
      WedgeTex = true;
    }
    if( (mask0 & tri::io::Mask::IOM_VERTTEXCOORD) && colorread) {
      m.vert.EnableTexCoord();
      VertTex = true;
    }
    if( (mask0 & tri::io::Mask::IOM_VERTNORMAL)) {
      //Rprintf("norm");
      //m.vert.EnableNormal();
      hasNormal = true;
    }
    if( (mask0 & tri::io::Mask::IOM_VERTQUALITY)) {
      m.vert.EnableQuality();
      hasQuality = true;
    }
    if( (mask0 & tri::io::Mask::IOM_FACEQUALITY)) {
      m.face.EnableQuality();
      hasFaceQuality = true;
    }
    tri::io::Mask::ClampMask(m, mask0);
    int err2 = 0;
    if (type == 0)
      err2 = tri::io::Importer<MyMesh>::Open(m,filename,mask0);
    else if (type == 1)
      err2 = tri::io::ImporterOFF<MyMesh>::Open(m,filename,mask0);
    if (err2) {
      return wrap(1);
    } else { 
      if (m.fn == 0) {
	updateNormals = false;
	clean = false;
      }
      
      
      // do texture processing
      bool tex = false;
      std::vector<float> texvec;
      std::vector<string> texfile;
      if (m.textures.size() > 0 && colorread) {

	if (!silent && clean)
	  Rprintf("To avoid wrong assignment of texture, cleaning has been disabled\n");
	
	clean = false;
	if (!VertTex && WedgeTex) {
	  m.vert.EnableTexCoord();	
	  tri::AttributeSeam::SplitVertex(m, ExtractVertex, CompareVertex);
	  vcg::tri::Allocator< MyMesh >::CompactVertexVector(m);
	  vcg::tri::Allocator< MyMesh >::CompactFaceVector(m);
	  texfile = m.textures;
	  texvec.resize(2*m.vn);
	}
	VertexIterator vi=m.vert.begin();
	for (int i=0;  i < m.vn; i++) {
	  texvec[i*2] = (*vi).T().U();
	  texvec[i*2+1] = (*vi).T().V();
	  vi++;
	}
      }
      
      if (clean) {
	int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
	int dupface = tri::Clean<MyMesh>::RemoveDuplicateFace(m);
	  int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
	vcg::tri::Allocator< MyMesh >::CompactVertexVector(m);
	vcg::tri::Allocator< MyMesh >::CompactFaceVector(m);
	if ((dup > 0 || unref > 0 || dupface > 0) && !silent)
	  Rprintf("Removed %i duplicate %i unreferenced vertices and %i duplicate faces\n",dup,unref,dupface);
      }  
      //setup indices
      SimpleTempData<MyMesh::VertContainer,int> indices(m.vert);
      // setup output structures
      NumericVector vb(3*m.vn);    
      std::vector<int> colvec;
      if (colorread && HasPerVertexColor(m)) {
	colvec.resize(3*m.vn);
      }
      IntegerVector it(3*m.fn);
      std::vector<double> normals;
      if (updateNormals || hasNormal) {
	normals.resize(3*m.vn);
      }
      std::vector<double> quality;
      if (hasQuality) {
	quality.resize(m.vn);
      }
      std::vector<double> facequality;
      if (hasFaceQuality) {
	facequality.resize(m.fn);
      }
    
      if (updateNormals) { // update Normals
	if (!hasNormal)
	  //m.vert.EnableNormal();
	tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
      }
      if (hasNormal && !updateNormals)
	tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
      // write back
      VertexIterator vi=m.vert.begin();
      for (int i=0;  i < m.vn; i++) {
	vb[i*3] = (*vi).P()[0];
	vb[i*3+1] = (*vi).P()[1];
	vb[i*3+2] = (*vi).P()[2];
	indices[vi] = i;
	if (updateNormals || hasNormal) {
	  normals[i*3] = (*vi).N()[0];
	  normals[i*3+1] = (*vi).N()[1];
	  normals[i*3+2] = (*vi).N()[2];
	}
	if (colorread && HasPerVertexColor(m)) {
	  colvec[i*3] = (*vi).C()[0];
	  colvec[i*3+1] = (*vi).C()[1];
	  colvec[i*3+2] = (*vi).C()[2];
	}
	if (hasQuality) {
	  quality[i] = (*vi).Q();
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
	    it[i*3]=vv[0];
	    it[i*3+1]=vv[1];
	    it[i*3+2]=vv[2];
	    if (hasFaceQuality) {
	        facequality[i] =(*fi).Q();
	    }
	    ++fi;
	  }
	}
      }
      
     
    
      return List::create(Named("vb") = vb, 
			  Named("it") = it,
			  Named("normals") = normals,
			  Named("colors") = colvec,
			  Named("texcoord") = texvec,
			  Named("texfile") = texfile,
			  Named("quality") = quality,
			  Named("facequality") = facequality
			  );
    }
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }

}
