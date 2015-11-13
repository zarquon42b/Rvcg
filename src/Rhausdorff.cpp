
#include "typedefHausdorff.h"


#include "RvcgIO.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace vcg;

std::string SaveFileName(const std::string &filename)
{
 int pos=filename.find_last_of('.',filename.length());
 std::string fileout=filename.substr(0,pos)+"_metro.ply";
 return fileout;
}

RcppExport SEXP Rhausdorff( SEXP vb0_, SEXP it0_,SEXP vb1_, SEXP it1_, SEXP vertSamp_, SEXP edgeSamp_, SEXP faceSamp_, SEXP unrefVert_, SEXP samplingType_, SEXP nSamples_, SEXP nSamplesArea_, SEXP saveMesh_, SEXP from_, SEXP to_, SEXP writeHist_, SEXP searchStruct_)
{
	try {
		// Read files
		CMesh m0, m1;
		Rvcg::IOMesh<CMesh>::RvcgReadR(m0,vb0_,it0_);
		Rvcg::IOMesh<CMesh>::RvcgReadR(m1,vb1_,it1_);
		// Declare variables
		bool vertSamp = Rcpp::as<bool>(vertSamp_);
		bool edgeSamp = Rcpp::as<bool>(edgeSamp_);
		bool faceSamp = Rcpp::as<bool>(faceSamp_);
		bool unrefVert = Rcpp::as<bool>(unrefVert_);
		unsigned int samplingType = Rcpp::as<unsigned int>(samplingType_);
		unsigned long nSamples = Rcpp::as<unsigned long>(nSamples_);
		double nSamplesArea = Rcpp::as<double>(nSamplesArea_);
		bool saveMesh = Rcpp::as<bool>(saveMesh_);
		double from = Rcpp::as<double>(from_);
		double to = Rcpp::as<double>(to_);
		bool writeHist = Rcpp::as<bool>(writeHist_);
		unsigned int searchStruct = Rcpp::as<unsigned int>(searchStruct_);

		float ColorMin=0, ColorMax=0;
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
			SamplingFlags::SIMILAR_SAMPLING;

		if (writeHist) flags |= SamplingFlags::HIST;
		if (vertSamp) flags &= ~SamplingFlags::VERTEX_SAMPLING;
		if (edgeSamp) flags &= ~SamplingFlags::EDGE_SAMPLING;
		if (faceSamp) flags &= ~SamplingFlags::FACE_SAMPLING;
		if (unrefVert) flags |= SamplingFlags::INCLUDE_UNREFERENCED_VERTICES;

		switch(samplingType){
			case 0 : flags = (flags | SamplingFlags::MONTECARLO_SAMPLING  ) & (~ SamplingFlags::NO_SAMPLING ); break;
			case 1 : flags = (flags | SamplingFlags::SUBDIVISION_SAMPLING ) & (~ SamplingFlags::NO_SAMPLING ); break;
			case 2 : flags = (flags | SamplingFlags::SIMILAR_SAMPLING     ) & (~ SamplingFlags::NO_SAMPLING ); break;
			default : Rprintf("%s\n","samplingType unknown" );
			exit(0);
		}

		if (nSamples != 0) NumberOfSamples = true; n_samples_target = nSamples;
		if (nSamplesArea != 0 ) SamplesPerAreaUnit = true; n_samples_per_area_unit = nSamplesArea;
		if (saveMesh) flags |= SamplingFlags::SAVE_ERROR;

		if (from != 0 && to != 0)
		{
			ColorMin=float(from); ColorMax=float(to);
		}

		switch(searchStruct){
			case 0 : flags |= SamplingFlags::USE_AABB_TREE; break;
			case 1 : flags |= SamplingFlags::USE_STATIC_GRID; break;
			case 2 : flags |= SamplingFlags::USE_HASH_GRID; break;
			case 3 : flags |= SamplingFlags::USE_OCTREE; break;
			default : Rprintf("%s\n","searchStruct unknown" );
			exit(0);
		}

		if(!(flags & SamplingFlags::USE_HASH_GRID) && !(flags & SamplingFlags::USE_AABB_TREE) && !(flags & SamplingFlags::USE_OCTREE)){
			flags |= SamplingFlags::USE_STATIC_GRID;
		}
		string m0NewName=SaveFileName("mesh1.ply");
		string m1NewName=SaveFileName("mesh2.ply");

		if(!NumberOfSamples && !SamplesPerAreaUnit){
				NumberOfSamples = true;
				n_samples_target = 10 * max(m0.fn,m1.fn);// take 10 samples per face
		}

		// compute face information
		tri::UpdateComponentEP<CMesh>::Set(m0);
		tri::UpdateComponentEP<CMesh>::Set(m1);

		// // set bounding boxes for S1 and S2
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
		unsigned long int nvertsamples0 = ForwardSampling.GetNVertexSamples();
		unsigned long int nsamples0 = ForwardSampling.GetNSamples();

		BackwardSampling.Hausdorff();
		dist1_max  = BackwardSampling.GetDistMax();
		double mean1 = BackwardSampling.GetDistMean();
		double rms1 = BackwardSampling.GetDistRMS();
		double nvertsamples1 = BackwardSampling.GetNVertexSamples();
		double nsamples1 = BackwardSampling.GetNSamples();

		// save error files.
		if(flags & SamplingFlags::SAVE_ERROR){
			vcg::tri::io::PlyInfo p;
			p.mask|=vcg::tri::io::Mask::IOM_VERTCOLOR | vcg::tri::io::Mask::IOM_VERTQUALITY /* | vcg::ply::PLYMask::PM_VERTQUALITY*/ ;
			//p.mask|=vcg::ply::PLYMask::PM_VERTCOLOR|vcg::ply::PLYMask::PM_VERTQUALITY;
			if(ColorMax!=0 || ColorMin != 0){
			vcg::tri::UpdateColor<CMesh>::PerVertexQualityRamp(m0,ColorMin,ColorMax);
			vcg::tri::UpdateColor<CMesh>::PerVertexQualityRamp(m1,ColorMin,ColorMax);
			}
			tri::io::ExporterPLY<CMesh>::Save( m0,m0NewName.c_str(),true,p);
			tri::io::ExporterPLY<CMesh>::Save( m1,m1NewName.c_str(),true,p);
		}

		// save error files.
		if(flags & SamplingFlags::HIST){
		  ForwardSampling.GetHist().FileWrite("forward_result.csv");
		  BackwardSampling.GetHist().FileWrite("backward_result.csv");
		}
		Rcpp::List out;
		out["ForwardSampling"] = Rcpp::List::create(
				Rcpp::Named("maxdist") = dist0_max,
				Rcpp::Named("meandist") = mean0,
				Rcpp::Named("RMSdist") = rms0,
				Rcpp::Named("nvbsamples") = nvertsamples0,
				Rcpp::Named("nsamples") = nsamples0
			);
		out["BackwardSampling"] = Rcpp::List::create(
				Rcpp::Named("maxdist") = dist1_max,
				Rcpp::Named("meandist") = mean1,
				Rcpp::Named("RMSdist") = rms1,
				Rcpp::Named("nvbsamples") = nvertsamples1,
				Rcpp::Named("nsamples") = nsamples1
			);
		return out;
	}
	catch (std::exception& e) {
		::Rf_error( e.what());
		return wrap(1);
	}
	catch (...) {
	::Rf_error("unknown exception");
	}
}