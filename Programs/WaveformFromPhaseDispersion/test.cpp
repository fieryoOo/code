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


int main() {
	std::string fphv("TEST/SyntheticFTAN/phv.txt");
	Curve<Point> phvs(fphv);
	SacRec sac;
	auto& shd = sac.shd;
	shd.delta = 0.0001;
	shd.b = 0.;
	shd.npts = (int)(0.5/shd.delta + 0.5);
	shd.e = shd.b + shd.delta * shd.npts;
	shd.dist = 500.;
	sac.sig.reset(new float[shd.npts]);
	float *sigsac = sac.sig.get();
	float twopi = M_PI * 2;
	double intpart;
	for(int i=0; i<shd.npts; i++) {
		float f = shd.b + shd.delta * i;
		float per = 1. / f;
		float phv = phvs.Val(per);
		float trvT = shd.dist / phv;
		float pha = modf(trvT / per, &intpart) - 0.125;
		if( pha > 0.5 ) pha -= 1.;
		pha *= twopi;
		sigsac[i] = pha;
	}
	sac.Write("test.SAC.ph");
}
