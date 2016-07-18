
#include "typedefMetro.h"


#include "RvcgIO.h"
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace vcg;
typedef  CMeshMetro::VertexIterator VertexIterator;


RcppExport SEXP Rmetro( SEXP mesh0_, SEXP mesh1_, SEXP vertSamp_, SEXP edgeSamp_, SEXP faceSamp_, SEXP unrefVert_, SEXP samplingType_, SEXP nSamples_, SEXP nSamplesArea_, SEXP from_, SEXP to_, SEXP searchStruct_, SEXP colormeshes_, SEXP silent_)
{
  try {
    CMeshMetro m0, m1;
    Rvcg::IOMesh<CMeshMetro>::mesh3d2Rvcg(m0,mesh0_);
    Rvcg::IOMesh<CMeshMetro>::mesh3d2Rvcg(m1,mesh1_);
    // Declare variables
    bool vertSamp = Rcpp::as<bool>(vertSamp_);
    bool edgeSamp = Rcpp::as<bool>(edgeSamp_);
    bool faceSamp = Rcpp::as<bool>(faceSamp_);
    bool unrefVert = Rcpp::as<bool>(unrefVert_);
    unsigned int samplingType = Rcpp::as<unsigned int>(samplingType_);
    unsigned long nSamples = Rcpp::as<unsigned long>(nSamples_);
    double nSamplesArea = Rcpp::as<double>(nSamplesArea_);
    double from = Rcpp::as<double>(from_);
    double to = Rcpp::as<double>(to_);
    unsigned int searchStruct = Rcpp::as<unsigned int>(searchStruct_);
    bool colormeshes = as<bool>(colormeshes_);
    bool silent = as<bool>(silent_);
    unsigned long n_samples_target;
    double n_samples_per_area_unit;

    // Sampling parameters
    double dist1_max, dist0_max;
    int flags;
    bool NumberOfSamples = false;
    bool SamplesPerAreaUnit = false;

    flags = SamplingFlags::VERTEX_SAMPLING |
      SamplingFlags::EDGE_SAMPLING |
      SamplingFlags::FACE_SAMPLING |
      //SamplingFlags::SIMILAR_SAMPLING |
      SamplingFlags::HIST;
      

    if (silent) flags |= SamplingFlags::SILENT;
    if (!vertSamp) flags &= ~SamplingFlags::VERTEX_SAMPLING;
    if (!edgeSamp) flags &= ~SamplingFlags::EDGE_SAMPLING;
    if (!faceSamp) flags &= ~SamplingFlags::FACE_SAMPLING;
    if (unrefVert) flags |= SamplingFlags::INCLUDE_UNREFERENCED_VERTICES;

    switch(samplingType){
    case 0 : flags |= SamplingFlags::MONTECARLO_SAMPLING  ; break;
    case 1 : flags |=  SamplingFlags::SUBDIVISION_SAMPLING; break;
    case 2 : flags |=  SamplingFlags::SIMILAR_SAMPLING ; break;
      //case 3 : flags |=  SamplingFlags::NO_SAMPLING ; break;
     default : ::Rf_error("%s\n","samplingType unknown" );
      
    }

    if (nSamples != 0) NumberOfSamples = true; n_samples_target = nSamples;
    if (nSamplesArea != 0 ) SamplesPerAreaUnit = true; n_samples_per_area_unit = nSamplesArea;
    flags |= SamplingFlags::SAVE_ERROR;

   
    switch(searchStruct){
    case 0 : flags |= SamplingFlags::USE_AABB_TREE; break;
    case 1 : flags |= SamplingFlags::USE_STATIC_GRID; break;
    case 2 : flags |= SamplingFlags::USE_HASH_GRID; break;
    case 3 : flags |= SamplingFlags::USE_OCTREE; break;
      //default : ::Rf_error("%s\n","searchStruct unknown" );
     
    }

    if(!(flags & SamplingFlags::USE_HASH_GRID) && !(flags & SamplingFlags::USE_AABB_TREE) && !(flags & SamplingFlags::USE_OCTREE)){
      flags |= SamplingFlags::USE_STATIC_GRID;
    }
		
    if(!NumberOfSamples && !SamplesPerAreaUnit){
      NumberOfSamples = true;
      n_samples_target = 10 * max(m0.fn,m1.fn);
      if (!silent && faceSamp)
	Rprintf("Number of samples set to %i\n",n_samples_target);// take 10 samples per face
    }

    // compute face information
    tri::UpdateComponentEP<CMeshMetro>::Set(m0);
    tri::UpdateComponentEP<CMeshMetro>::Set(m1);

    // // set bounding boxes for S1 and S2
    tri::UpdateBounding<CMeshMetro>::Box(m0);
    tri::UpdateBounding<CMeshMetro>::Box(m1);

    // set Bounding Box.
    Box3<CMeshMetro::ScalarType>    bbox, tmp_bbox_M1=m0.bbox, tmp_bbox_M2=m1.bbox;
    bbox.Add(m0.bbox);
    bbox.Add(m1.bbox);
    bbox.Offset(bbox.Diag()*0.02);
    m0.bbox = bbox;
    m1.bbox = bbox;

    // sampling
    Sampling<CMeshMetro> ForwardSampling(m0,m1);
    Sampling<CMeshMetro> BackwardSampling(m1,m0);

    ForwardSampling.SetFlags(flags);
    if(NumberOfSamples){
      ForwardSampling.SetSamplesTarget(n_samples_target);
      n_samples_per_area_unit = ForwardSampling.GetNSamplesPerAreaUnit();
    }
    else{
      ForwardSampling.SetSamplesPerAreaUnit(n_samples_per_area_unit);
      n_samples_target = ForwardSampling.GetNSamplesTarget();
    }

    BackwardSampling.SetFlags(flags);
    if(NumberOfSamples)
      {
	BackwardSampling.SetSamplesTarget(n_samples_target);
	n_samples_per_area_unit = BackwardSampling.GetNSamplesPerAreaUnit();
      }
    else{
      BackwardSampling.SetSamplesPerAreaUnit(n_samples_per_area_unit);
      n_samples_target = BackwardSampling.GetNSamplesTarget();
    }

    ForwardSampling.Hausdorff();
    dist0_max  = ForwardSampling.GetDistMax();
    double mean0 = ForwardSampling.GetDistMean();
    double rms0 = ForwardSampling.GetDistRMS();
    double area0 = ForwardSampling.GetArea();
    double distvol0 = ForwardSampling.GetDistVolume();
    unsigned long int nvertsamples0 = ForwardSampling.GetNVertexSamples();
    unsigned long int nsamples0 = ForwardSampling.GetNSamples();

    BackwardSampling.Hausdorff();
    dist1_max  = BackwardSampling.GetDistMax();
    double mean1 = BackwardSampling.GetDistMean();
    double rms1 = BackwardSampling.GetDistRMS();
    double area1 = BackwardSampling.GetArea();
    double distvol1 = BackwardSampling.GetDistVolume();
    double nvertsamples1 = BackwardSampling.GetNVertexSamples();
    double nsamples1 = BackwardSampling.GetNSamples();

    //write heatmap to vertex color
    if(colormeshes){
      if(from != 0 || to != 0){
	vcg::tri::UpdateColor<CMeshMetro>::PerVertexQualityRamp(m0,from,to);
	vcg::tri::UpdateColor<CMeshMetro>::PerVertexQualityRamp(m1,from,to);
      } else {
	vcg::tri::UpdateColor<CMeshMetro>::PerVertexQualityRamp(m0);
	vcg::tri::UpdateColor<CMeshMetro>::PerVertexQualityRamp(m1);
      }
    }
    // write back color information
    std::vector<int> colvec0(3*m0.vn),colvec1(3*m1.vn);
    std::vector<float> quality0(m0.vn), quality1(m1.vn);
    
    VertexIterator vi=m0.vert.begin();
    for (int i=0;  i < m0.vn; i++) {
      colvec0[i*3] = (*vi).C()[0];
      colvec0[i*3+1] = (*vi).C()[1];
      colvec0[i*3+2] = (*vi).C()[2];
      quality0[i] = (*vi).Q();
      ++vi;
    }
    vi=m1.vert.begin();
    for (int i=0;  i < m1.vn; i++) {
      colvec1[i*3] = (*vi).C()[0];
      colvec1[i*3+1] = (*vi).C()[1];
      colvec1[i*3+2] = (*vi).C()[2];
      quality1[i] = (*vi).Q();
      ++vi;
    }
    
    List mesh0 = Rvcg::IOMesh<CMeshMetro>::RvcgToR(m0);
    //mesh0["quality"] = quality0;
    List mesh1 = Rvcg::IOMesh<CMeshMetro>::RvcgToR(m1);
    //mesh1["quality"] = quality1;
    // save error files.
		
    // create sampling histogram as R matrices
    Histogram<double> fwdhist = ForwardSampling.GetHist();
    int nfwd = fwdhist.BinNum();
    double fwdcnt = fwdhist.Cnt();
    NumericMatrix forward_hist(nfwd+2,2);
    for (int i = 0; i <= nfwd+1; i++) {
      double lbi = fwdhist.BinLowerBound(i);
      double hi = fwdhist.BinCountInd(i)/fwdcnt;
      forward_hist(i,0) = fwdhist.BinLowerBound(i);
      forward_hist(i,1) = fwdhist.BinCountInd(i)/fwdcnt;
    }
    Histogram<double> bckhist = BackwardSampling.GetHist();
    int nbck = bckhist.BinNum();
    double bckcnt = bckhist.Cnt();
    NumericMatrix backward_hist(nbck+2,2);
    for (int i = 0; i <= nbck+1; i++) {
      double lbi = bckhist.BinLowerBound(i);
      double hi = bckhist.BinCountInd(i)/bckcnt;
      backward_hist(i,0) = bckhist.BinLowerBound(i);
      backward_hist(i,1) = bckhist.BinCountInd(i)/bckcnt;
    }
		    
		  
		
    Rcpp::List out;
    out["ForwardSampling"] = Rcpp::List::create(
						Rcpp::Named("maxdist") = dist0_max,
						Rcpp::Named("meandist") = mean0,
						Rcpp::Named("RMSdist") = rms0,
						Rcpp::Named("area") = area0,
						Rcpp::Named("distvolume") = distvol0,
						Rcpp::Named("nvbsamples") = nvertsamples0,
						Rcpp::Named("nsamples") = nsamples0
						);
    out["BackwardSampling"] = Rcpp::List::create(
						 Rcpp::Named("maxdist") = dist1_max,
						 Rcpp::Named("meandist") = mean1,
						 Rcpp::Named("RMSdist") = rms1,
						 Rcpp::Named("area") = area1,
						 Rcpp::Named("distvolume") = distvol1,
						 Rcpp::Named("nvbsamples") = nvertsamples1,
						 Rcpp::Named("nsamples") = nsamples1
						 );
    out["forward_hist"] = forward_hist;
    out["backward_hist"] = backward_hist;
    out["distances1"] = quality0;
    out["distances2"] = quality1;
    if (colormeshes) {
    out["mesh1"] = Rcpp::List::create(
				      Rcpp::Named("mesh") = mesh0,
				      Rcpp::Named("colors") = colvec0
				      );
    out["mesh2"] = Rcpp::List::create(
				      Rcpp::Named("mesh") = mesh1,
				      Rcpp::Named("colors") = colvec1
				      );
    }
    
    
    return out;
  }
  catch (std::exception& e) {
    ::Rf_error( e.what());
  }
  catch (...) {
    ::Rf_error("unknown exception");
  }
}
