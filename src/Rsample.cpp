#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
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


RcppExport SEXP Rsample(SEXP mesh_, SEXP SampleNum_, SEXP type_, SEXP MCsamp_, SEXP geodes_)
{
  try {// declare Mesh and helper variables
    int SampleNum = Rcpp::as<int>(SampleNum_);  
    //double tol = Rcpp::as<double>(tol_);  
    const int type = Rcpp::as<int>(type_);  
    const int MCsamp = Rcpp::as<int>(MCsamp_);
    const bool geodes = Rcpp::as<bool>(geodes_);
    
    MyMesh m,msamp;
    float radius;
    VertexIterator vi;
    FaceIterator fi;
    // allocate mesh and fill it
    int check = Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(m,mesh_);
    if (check != 0) {
      return wrap(1);
    } else {
      vector<Point3f> myVec;
      typedef TrivialSampler<MyMesh>  BaseSampler ;
 
      BaseSampler ts2(myVec);
      std::vector<Point3f> MontecarloSamples, poissonsamples;
      BaseSampler mcSampler(MontecarloSamples);
      BaseSampler pdSampler(poissonsamples);
  
  
      if (type == 1) {
	tri::MontecarloSampling(m,poissonsamples,SampleNum);
      } else { //this is poison disk sampling with more options than the wrapper commented below
	radius = tri::SurfaceSampling<MyMesh,BaseSampler>::ComputePoissonDiskRadius(m,SampleNum);
	tri::SurfaceSampling<MyMesh, BaseSampler>::PoissonDiskParam pp;
	pp.preGenMesh=&m;
	//tri::SurfaceSampling<MyMesh, BaseSampler>::PoissonDiskParam::Stat stat;
	//pp.pds = &stat;
	//pp.pds->sampleNum = SampleNum;
	pp.geodesicDistanceFlag=geodes;
	MyMesh MontecarloMesh;
	SurfaceSampling<MyMesh,BaseSampler>::Montecarlo(m, mcSampler, SampleNum*MCsamp);
    
	tri::Allocator<MyMesh>::AddVertices(MontecarloMesh,MontecarloSamples.size());
	for(unsigned int j=0;j < MontecarloSamples.size();++j)
	  MontecarloMesh.vert[j].P()=MontecarloSamples[j];
    
	tri::UpdateBounding<MyMesh>::Box(MontecarloMesh);
	tri::SurfaceSampling<MyMesh,BaseSampler>::PoissonDiskPruning(pdSampler, MontecarloMesh, radius, pp);
      }
      unsigned int outsize = poissonsamples.size();
      Rcpp::NumericMatrix vbout(3,outsize);
  
      for (unsigned int i=0;  i < outsize; i++) {
	Point3f tmp = poissonsamples[i];
	vbout(0,i) = tmp[0];
	vbout(1,i) =  tmp[1];
	vbout(2,i) =  tmp[2];
      }
  
      return Rcpp::wrap(vbout);
    }
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
 

    
