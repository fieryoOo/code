#include "SacRec.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

struct Ddata {

	Ddata( float perin=NaN, float velin=NaN ) 
		: per(perin), vel(velin) {}

	Ddata(const std::string& input) {
		int nrd = sscanf(input.c_str(), "%f %f", &per, &vel);
		if( nrd < 2 )
			throw std::runtime_error("Error(Ddata::Ddata): format error");
	}

	friend bool operator< ( const Ddata& d1, const Ddata& d2 ) {
		return (d1.per < d2.per);
	}

	friend std::ostream& operator<< ( std::ostream& o, const Ddata& ddata ) {
		o << ddata.per << " " << ddata.vel;
		return o;
	}

	static constexpr float NaN = -12345.;
	float per = NaN, vel = NaN;
};


class DISP {
public:
		DISP( const std::string& fname ) { Load(fname); }
		void Load( const std::string& fname ) {
			// check infile
			std::ifstream fin(fname);
			if( ! fin )
				throw std::runtime_error("Error(DISP::Load): cannot access dispersion infile " + fname);
			// read
			for(std::string line; std::getline(fin, line); ) {
				dataV.push_back( Ddata(line) );
			}
			std::cout<<dataV.size()<<" points loaded from file "<<fname<<std::endl;
			// sort
			std::sort(dataV.begin(), dataV.end());
			// debug
			//for( const auto& d : dataV ) std::cerr<<d<<std::endl;
		}

		float permin() { return dataV.front().per; }
		float permax() { return dataV.back().per; }

		float Vel( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), Ddata(per) );
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return Ddata::NaN;
			const auto itlow = itupp - 1;
			float vdiff = (*itupp).vel - (*itlow).vel, Tdiff = (*itupp).per - (*itlow).per;
			float deriv = vdiff / Tdiff;
			return (*itlow).vel + Tdiff * deriv;
		}

		float Derivative( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), Ddata(per) );
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return Ddata::NaN;
			const auto itlow = itupp - 1;
			return ( (*itupp).vel - (*itlow).vel ) / ( (*itupp).per - (*itlow).per );
		}

private:
		std::vector<Ddata> dataV;
};

int main(int argc, char* argv[]) {
	if( argc != 4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [phv infile] [grv outfile] [sacout]"<<std::endl;
		exit(-1);
	}
	
	// read phase dispersion
	DISP disp(argv[1]);

	// synthetic at center period 10 sec and distance 100km
	int npts = 10000;
	float delta = 1.;
	float pi=3.14159265, dom = 2.*pi*0.1, oopi = 1./pi;
	float dist = 500.;
	SacRec sacout;
	sacout.sig.reset( new float[npts]() );
	sacout.shd.b = 0.;
	sacout.shd.dist = dist;
	sacout.shd.delta = delta;
	sacout.shd.npts = npts;
	std::vector<Ddata> grvV;
	std::cout<<" per range: "<<disp.permin()<<" - "<<disp.permax()<<std::endl;
	float perfactor = exp(0.005);	// log(per) + 0.005
	for( float per=disp.permin(); per<disp.permax(); per*=perfactor ) {
	//for( float per=5.; per<10.; per*=perfactor ) {
		float c = disp.Vel(per), om0 = 2.*pi/per;
		float wavenum = 2.*pi / (per * c);
		float kderiv = 1./c + per/(c*c)*disp.Derivative(per);
		grvV.push_back( Ddata(per, 1/kderiv) );
		float* sigsac = sacout.sig.get();
		for(int it=0; it<npts; it++) {
			float t = it * delta;
			float Y = 0.5*dom * (t - kderiv*dist);
			sigsac[it] += dom*oopi * sin(Y)/Y * cos(om0*t - wavenum*dist);
		}
	}

	// output group dispersion
	std::ofstream fout(argv[2]);
	if( ! fout ) {
		std::cerr<<"Error(main): cannot write to file "<<argv[2]<<std::endl;;
		exit(-2);
	}
	for(int iper=0; iper<grvV.size(); iper++)
		fout<<grvV[iper]<<std::endl;

	// output sac
	sacout.Write( argv[3] );
	

	return 0;
}
