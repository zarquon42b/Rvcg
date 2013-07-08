#include <../typedef.h>
#include <../RvcgIO.h>
#include <Rcpp.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <limits>
#include <vcg/complex/algorithms/clustering.h>
#include <vcg/simplex/face/distance.h>
#include <vcg/complex/algorithms/geodesic.h>
using namespace vcg;
using namespace tri;
using namespace Rcpp;
using namespace std;


RcppExport SEXP Rsample(SEXP _vb, SEXP _it, SEXP _SampleNum, SEXP _type)
{
  // declare Mesh and helper variables
  int SampleNum = Rcpp::as<int>(_SampleNum);  
  //double tol = Rcpp::as<double>(_tol);  
  const int type = Rcpp::as<int>(_type);  
  int i, j;
  MyMesh m,msamp;
  float radius = 0;
  VertexIterator vi;
  FaceIterator fi;
  // allocate mesh and fill it
  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,_vb,_it);
  m.vert.EnableVFAdjacency();
  m.face.EnableFFAdjacency();
  m.face.EnableVFAdjacency();
  vector<Point3f> myVec;
  typedef TrivialSampler<MyMesh>  BaseSampler ;
  BaseSampler mcSample(myVec);
   BaseSampler ts2(myVec);

  
  
  //SurfaceSampling<MyMesh, TrivialSampler<MyMesh> >::Poissondi(m, ts);
  Rcpp::NumericMatrix vbout(3,SampleNum);
  /*
  MyMesh MontecarloMesh;
  vcg::tri::Append<MyMesh,MyMesh>::Mesh(MontecarloMesh,m);
  MyMesh *presampledmesh;
  tri::SurfaceSampling<MyMesh,BaseSampler>::Montecarlo(MontecarloMesh, mcSample, SampleNum*5);
  */
  std::vector<vcg::Point3f> ExactVec;
  std::vector<vcg::Point3f> PerturbVec;
  if (type == 1)
    tri::MontecarloSampling(m,ExactVec,SampleNum);
  else
    vcg::tri::PoissonSampling(m, ExactVec, SampleNum,radius);
    //tri::SurfaceSampling<MyMesh,BaseSampler>::PoissonDiskPruning(ts2, m, radius);
  //tri::PoissonDisk(m,ExactVec,10);
  for (i=0;  i < SampleNum; i++) 
    {
      Point3f tmp = ExactVec[i];
      vbout(0,i) = tmp[0];
      vbout(1,i) =  tmp[1];
      vbout(2,i) =  tmp[2];
    }
	 
	  return Rcpp::wrap(vbout);
	  /*return Rcpp::List::create(Rcpp::Named("vb") = vbout,
			    Rcpp::Named("it") = itout,
			    Rcpp::Named("normals") = normals
			    );
	  */
}
 

    
