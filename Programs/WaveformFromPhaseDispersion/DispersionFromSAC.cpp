#include "SacRec.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

/* ----- single point data structures ----- */
template <class T> class Curve;
class KDeriv;
class Point {
public:
	Point( float perin=NaN, float velin=NaN, float sdensin=1. ) 
		: per(perin), vel(velin), sdensity(sdensin) {}

	Point( const std::string& input ) {
		int nrd = sscanf(input.c_str(), "%f %f %f", &per, &vel, &sdensity);
		if( nrd < 2 )
			throw std::runtime_error("Error(Point::Point): format error");
	}

	friend bool operator< ( const Point& p1, const Point& p2 ) {
		return (p1.per < p2.per);
	}

	friend std::ostream& operator<< ( std::ostream& o, const Point& pt ) {
		o << pt.per << " " << pt.vel << " " << pt.sdensity;
		return o;
	}

	friend Curve<Point>;
	friend Curve<Point> operator-(const Curve<Point>& c1, const Curve<Point>& c2);

	static constexpr float NaN = -12345., pi = 3.1415926536;

protected:
	float per = NaN, vel = NaN, sdensity = 1.;
};

/* ----- data containers ----- */
template <class T>
class Curve {
public:
		Curve() {}
		Curve( const std::string& fname ) { Load(fname); }
		void Load( const std::string& fname ) {
			// check infile
			std::ifstream fin(fname);
			if( ! fin )
				throw std::runtime_error("Error(Curve::Load): cannot access dispersion infile " + fname);
			// read
			for(std::string line; std::getline(fin, line); ) {
				dataV.push_back( T(line) );
			}
			std::cout<<dataV.size()<<" points loaded from file "<<fname<<std::endl;
			// sort
			std::sort(dataV.begin(), dataV.end());
			// debug
			//for( const auto& d : dataV ) std::cerr<<d<<std::endl;
		}

		void clear() { dataV.clear();	}
		void push_back( float per, float vel, float om ){ dataV.push_back( Point(per, vel, om) ); }

		float xmin() { return dataV.front().per; }
		float xmax() { return dataV.back().per; }

		float Val( const float per ) const {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), T(per) );
			if( itupp == dataV.end() ) {
				if( fabs((*(itupp-1)).per - per) < 0.1 ) return (*(itupp-1)).vel;
				else return T::NaN;
			}
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return T::NaN;
			const auto itlow = itupp - 1;
			float vdiff = (*itupp).vel - (*itlow).vel, Tdiff = (*itupp).per - (*itlow).per;
			float deriv = vdiff / Tdiff;
			return (*itlow).vel + ( per-(*itlow).per ) * deriv;
		}

		float Deriv( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), T(per) );
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return T::NaN;
			const auto itlow = itupp - 1;
			return ( (*itupp).vel - (*itlow).vel ) / ( (*itupp).per - (*itlow).per );
		}

		void FindValue(float val, std::vector<float>& perV) {
			perV.clear();
			if( dataV.size() <= 1 ) return;
			for( int i=0; i<dataV.size()-1; i++ ) {
				const auto& datacur = dataV[i];
				float velh = dataV[i+1].vel, vell = datacur.vel;
				if( val == vell ) {	// equal to or
					perV.push_back(datacur.per);
				} else if( (val-vell) * (val-velh) < 0. ) { // in between
					float perh = dataV[i+1].per, perl = datacur.per;
					float vdiff = velh - vell, pdiff = perh - perl;
					//std::cerr<<"looking for "<<val<<" in between "<<perl<<"-"<<perh<<"sec "<<vell<<"-"<<velh<<"s/km"<<std::endl;
					perV.push_back( perl + (val-vell) * pdiff / vdiff );
				}
			}
			const auto& datacur = dataV.back();
			if( val == datacur.vel ) perV.push_back(datacur.per);
		}

		void Write( const std::string& fname ) {
			std::ofstream fout( fname );
			if( ! fout ) {
				std::cerr<<"Error(Curve::Write): cannot write to file "<<fname<<std::endl;;
				exit(-2);
			}
			for( const auto& spP : dataV )
				fout<<spP<<"\n";
		}

		friend Curve<Point> operator-(const Curve<Point>& c1, const Curve<Point>& c2) {
			Curve<Point> c3;
			for( auto p1 : c1.dataV ) {
				float val2 = c2.Val(p1.per);
				if( val2 == Point::NaN ) continue;
				p1.vel -= val2;
				c3.dataV.push_back(p1);
			}
			return c3;
		}

private:
		std::vector<T> dataV;
};

inline int nint( float f ) { return (int)floor(f+0.5); }
int main(int argc, char* argv[]) {
	if( argc != 3 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sac_in] [tmax]\n";
		return -1;
	}

	try {
		// load sac and suppress noise
		SacRec sac(argv[1]);
		sac.Load();
		//sac.cut(0., 1000.);
		float *sigsac = sac.sig.get();
		int tcbeg = nint( atof(argv[2]) ), tclen = 100;
		int ib = nint(tcbeg/sac.shd.delta), ie = nint((tcbeg+tclen)/sac.shd.delta), i = ib;
		float ftmp_cos = M_PI * 0.5 / (ie-ib);
		for( ; i<ie; i++ ) sigsac[i] *= cos( ftmp_cos * (i-ib) );
		for( ; i<sac.shd.npts; i++ ) sigsac[i] = 0.;

		// FFT
		SacRec sac_am, sac_ph;
		sac.ToAmPh( sac_am, sac_ph );
		float dist = sac.shd.dist;
		sac_am.clear(); sac.clear();

		// correct 2pi in phase
	sac_ph.Write("1.SAC");
		float twopi = 2. * M_PI;
		auto& shd = sac_ph.shd;
		sigsac = sac_ph.sig.get();
		for(int i=1; i<shd.npts; i++) {
			sigsac[i] = sigsac[i] - nint( (sigsac[i]-sigsac[i-1])/twopi ) * twopi;
		}
		SacRec sac_tmp;
		sac_ph.Smooth(0.01, sac_tmp);
		sac_ph = std::move(sac_tmp);
	sac_ph.Write("2.SAC");

		// output file
		std::string outname(argv[1]); outname += "_phv";
		std::ofstream fout( outname );
		if( ! fout ) {
			std::cerr<<"Error(main): cannot write to file "<<outname<<"\n";
			return -2;
		}

		// work on the phase
		sigsac = sac_ph.sig.get();
		float vel = 4.5;//, twopi = M_PI * 2.;
		int ibeg = nint(0.03/shd.delta);
		for(int i=ibeg; i<shd.npts; i++) {
			float freq = i*shd.delta;
			float ftmp1 = dist*freq, ftmp2 = sigsac[i]/twopi + 0.125;
			float n2pi = nint( ftmp1/vel - ftmp2 );
			//std::cerr<<ftmp1<<" "<<vel<<" "<<ftmp2<<" "<<n2pi<<std::endl;
			vel = ftmp1 / (ftmp2 + n2pi);
			fout<< 1./freq << " " << vel << "\n";
		}
	} catch(...) {
		std::cerr<<"Error(main): exception detected!\n";
	}

	return 0;
}
