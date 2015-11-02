
#include "typedefHausdorff.h"


#include "RvcgIO.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace vcg;


RcppExport SEXP Rhausdorff( SEXP vb0_, SEXP it0_,SEXP vb1_, SEXP it1_)
{
	try {
		// Read files
		CMesh m0, m1;
		Rvcg::IOMesh<CMesh>::RvcgReadR(m0,vb0_,it0_);
		Rvcg::IOMesh<CMesh>::RvcgReadR(m1,vb1_,it1_);

		// Sampling parameters
		double dist1_max, dist0_max;
		int flags;

		flags = SamplingFlags::VERTEX_SAMPLING |
		SamplingFlags::EDGE_SAMPLING |
		SamplingFlags::FACE_SAMPLING |
		SamplingFlags::SIMILAR_SAMPLING;

		// compute face information
		tri::UpdateComponentEP<CMesh>::Set(m0);
		tri::UpdateComponentEP<CMesh>::Set(m1);

		// set bounding boxes for S1 and S2
		tri::UpdateBounding<CMesh>::Box(m0);
		tri::UpdateBounding<CMesh>::Box(m1);

		// set Bounding Box.
		Box3<CMesh::ScalarType>    bbox, tmp_bbox_M1=m0.bbox, tmp_bbox_M2=m1.bbox;
		bbox.Add(m0.bbox);
		bbox.Add(m1.bbox);
		bbox.Offset(bbox.Diag()*0.02);
		m0.bbox = bbox;
		m1.bbox = bbox;

		// sampling
		Sampling<CMesh> ForwardSampling(m0,m1);
		Sampling<CMesh> BackwardSampling(m1,m0);

		double n_samples_target (100000);
		double SamplesPerAreaUnit;
		ForwardSampling.SetSamplesTarget(n_samples_target);
		SamplesPerAreaUnit = ForwardSampling.GetNSamplesPerAreaUnit();

		ForwardSampling.SetFlags(flags);
			// ForwardSampling.SetSamplesTarget(n_samples_target);
			// SamplesPerAreaUnit = ForwardSampling.GetNSamplesPerAreaUnit();
		ForwardSampling.SetSamplesPerAreaUnit(SamplesPerAreaUnit);
		ForwardSampling.Hausdorff();
		dist0_max  = ForwardSampling.GetDistMax();
		double mean0 = ForwardSampling.GetDistMean();
		double rms0 = ForwardSampling.GetDistRMS();
		double nvertsamples0 = ForwardSampling.GetNVertexSamples();
		double nsamples0 = ForwardSampling.GetNSamples();


		BackwardSampling.SetFlags(flags);
			// BackwardSampling.SetSamplesTarget(n_samples_target);
			// SamplesPerAreaUnit = BackwardSampling.GetNSamplesPerAreaUnit();
		BackwardSampling.SetSamplesPerAreaUnit(SamplesPerAreaUnit);
		BackwardSampling.Hausdorff();
		dist1_max  = BackwardSampling.GetDistMax();
		double mean1 = BackwardSampling.GetDistMean();
		double rms1 = BackwardSampling.GetDistRMS();
		double nvertsamples1 = BackwardSampling.GetNVertexSamples();
		double nsamples1 = BackwardSampling.GetNSamples();

		return Rcpp::List::create(
				Rcpp::Named("maxdist") = dist0_max,
				Rcpp::Named("meandist") = mean0,
				Rcpp::Named("RMSdist") = rms0,
				Rcpp::Named("nvbsamples") = nvertsamples0,
				Rcpp::Named("nsamples") = bbox.DimX()
			);

		// mesh1 = Rcpp::List::create(
		// 		Rcpp::Named("maxdist") = dist1_max,
		// 		Rcpp::Named("meandist") = mean1,
		// 		Rcpp::Named("RMSdist") = rms1,
		// 		Rcpp::Named("nvbsamples") = nvertsamples1,
		// 		Rcpp::Named("nsamples") = nsamples1,
		// 	);
		// return Rcpp::List::create(
		// 	Rcpp::Named("forward") = mesh0,
		// 	Rcpp::Named("backward") = mesh1,
		// 	);
	}
	catch (std::exception& e) {
		::Rf_error( e.what());
		return wrap(1);
	}
	catch (...) {
	::Rf_error("unknown exception");
	}
}