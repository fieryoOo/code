#include "Curve.h"
#include "VectorOperations.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>

class ModelFile {
public:
	const float Tsed() const { return thick[0]; }
	const float Tcru() const { return thick[1]; }
	const float Tman() const { return thick[2]; }
	

public:

	ModelFile( const std::string& fname = "", const size_t col_no = 2 ) {
		if( ! fname.empty() ) Load(fname, col_no);
	}

	// load data (colx=1, coly=col_no) from file
	void Load( const std::string& fname, const size_t col_no ) {
		std::ifstream fin( fname );
		if( ! fin )
			throw std::runtime_error("IO failed on " + fname);
		for( std::string line; std::getline(fin, line); ) {
			std::stringstream ss(line);
			PointC p;
			// x data
			bool suc = (ss >> p.x);
			// y data at col_no
			for( int icol=2; icol<=col_no; icol++ ) suc = (ss >> p.y);
			if( ! suc ) 
				throw std::runtime_error("requested column# " + std::to_string(col_no) + "out of range" );
			// Q at the last col
			float Q = p.y;
			while( !ss.eof() ) ss >> Q;
			// determine which layer the current slice belongs to
			int ilayer;
			if( Q == Qsed ) {
				ilayer = 0;
			} else if ( Q == Qcru ) {
				ilayer = 1;
			} else { 
				ilayer = 2; 
			}
			// compare current and the last layer
			if( ! (ilayer==ilayerlast || ilayer==ilayerlast+1) )
				throw std::runtime_error("unexpected layers/Q values");
			// push back
			dataV(ilayer).push_back(p);
			thick[ilayer] += p.x;
			ilayerlast=ilayer;
		}
		std::cout<<"### ModelFile::Load: Tsed="<<Tsed()<<"("<<npsed()<<") Tcru="<<Tcru()<<"("<<npcru()<<") Tman="<<Tman()<<"("<<npman()<<")"<<std::endl;
	}

	inline int npsed() const { return dataVsed.size(); }
	inline int npcru() const { return dataVcru.size(); }
	inline int npman() const { return dataVman.size(); }

	// normalize layer thicknesses by shrinking/expanding each model file
	void NormalizeThick(const float Tsed_avg, const float Tcru_avg, const float Tman_avg) {
		float nfactor[3] = { Tsed_avg/thick[0], Tcru_avg/thick[1], Tman_avg/thick[2] };
		for( int ilayer=0; ilayer<3; ilayer++ ) {
			auto& DataV = dataV(ilayer);
			// skip empty layer
			if( DataV.size() == 0 ) continue;
			// check layer thickness
			float normfactor = nfactor[ilayer];
			if( normfactor != normfactor ) 
				throw std::runtime_error("layer #"+std::to_string(ilayer)+" has 0 thickness");
			// normalize each slices
			for( auto& p : DataV ) p.x *= normfactor;
		}
		// update thickness
		thick[0] = Tsed_avg; thick[1] = Tcru_avg; thick[2] = Tman_avg;
	}

	Curve<PointC> ToCurve() const {
		Curve<PointC> c1; c1.FromThickData( dataVsed, 0., true, false );
		Curve<PointC> c2; c2.FromThickData( dataVcru, Tsed(), false, false );
		Curve<PointC> c3; c3.FromThickData( dataVman, Tsed()+Tcru(), false, true );
		c1.append(c2); c1.append(c3);
		return c1;
	}

protected:
	static constexpr float Qsed = 80.00;
	static constexpr float Qcru = 599.99;

private:
	int ilayerlast = 0;
	std::vector<PointC> dataVsed;
	std::vector<PointC> dataVcru;
	std::vector<PointC> dataVman;
	std::vector<PointC>& dataV( const int ilayer ) {
		switch(ilayer) {
			case 0: return dataVsed;
			case 1: return dataVcru;
			case 2: return dataVman;
			default: throw std::runtime_error("ilayer out of range");
		}
	}

	float thick[3] = {0., 0., 0.};
};




int main( int argc, char* argv[] ) {
	if( argc != 3 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [model file list] [model out file]"<<std::endl;
		return -1;
	}

	std::ifstream fin( argv[1] );
	if( ! fin ) {
		std::cerr<<"IO failed on "<<argv[1]<<std::endl;
		return -2;
	}

	// load all model files and compute average thicknesses
	float Tsed_avg = 0., Tcru_avg = 0., Tman_avg = 0.;
	std::vector<ModelFile> MFV;
	for( std::string line; std::getline(fin, line); ) {
		// load in thickness(col 1) and Vs(col 3)
		ModelFile mf(line, 3);
		Tsed_avg += mf.Tsed();
		Tcru_avg += mf.Tcru();
		Tman_avg += mf.Tman();
		MFV.push_back( mf );
	}
	const size_t Nmod = MFV.size();
	Tsed_avg /= Nmod; Tcru_avg /= Nmod; Tman_avg /= Nmod;
	std::cout<<"Average thicknesses: "<<Tsed_avg<<" "<<Tcru_avg<<" "<<Tman_avg<<std::endl;

	// normalize layer thicknesses by shrinking/expanding each model file
	for( auto& mf : MFV )
		mf.NormalizeThick(Tsed_avg, Tcru_avg, Tman_avg);
	// find the layer with the largest Sed_npts to start with
	int nmax = 0, imfmax = 0;
	for( int i=0; i<MFV.size(); i++ ) {
		const auto& mf = MFV[i];
		int npsed = mf.npsed();
		if( nmax<npsed ) { imfmax = i; nmax = npsed;	}
	}
	// store curves into vector
	std::vector< Curve<PointC> > cvV; cvV.reserve( MFV.size() );
	cvV.push_back( MFV[imfmax].ToCurve() );
	for( int i=0; i<MFV.size(); i++ ) {
		if( i == imfmax ) continue;
		cvV.push_back( MFV[i].ToCurve() );
	}
	// compute mean and std
	Curve<PointC> CVmean, CVstd;
	VO::MeanSTD( cvV.begin(), cvV.end(), CVmean, CVstd );
	// output results
	//CVmean.Write(argv[2]);
	//CVstd.Write(std::string(argv[2])+"_std");
	std::ofstream fout(argv[2]);
	PointC mean, std;
	for( CVmean.rewind(), CVstd.rewind(); CVmean.get(mean)&&CVstd.get(std); CVmean.next(), CVstd.next() )
		fout<<mean.x<<" "<<mean.y<<" "<<std.y<<"\n";

	return 0;
}

