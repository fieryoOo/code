#include "SacRec.h"
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

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
	friend KDeriv;

	static constexpr float NaN = -12345., pi = 3.1415926536;

protected:
	float per = NaN, vel = NaN, sdensity = 1.;
};

class Dispersion;
class Ddata : public Point {
public:
	Ddata( float perin=NaN, float velin=NaN, float sdensin=1. )
		: Point(perin, velin, sdensin) {}

	Ddata(const std::string& input) 
		: Point(input) {
		om = 2. * pi / per;
		k_ang = om / vel;
	}

	friend Curve<Ddata>;
	friend Dispersion;

private:
	float om, k_ang;	//angular wavenumber
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

		float Val( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), T(per) );
			if( itupp == dataV.end() ) {
				if( fabs((*(itupp-1)).per - per) < 0.1 ) return (*(itupp-1)).vel;
				else return Ddata::NaN;
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

protected:
		std::vector<T> dataV;
};

class KDeriv : public Curve<Point> {
public:
		float Deriv_2om( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), Point(per) );
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return Point::NaN;
			const auto itlow = itupp - 1;
			return ( (*itupp).vel - (*itlow).vel ) / ( (*itupp).sdensity - (*itlow).sdensity );
		}

		float Reciprocal( std::vector<Point>& grvV ) {
			grvV.resize( dataV.size() );
			for(int i=0; i<dataV.size(); i++)
				grvV[i] = Point(dataV[i].per, 1./dataV[i].vel, dataV[i].sdensity);
		}
};

class Dispersion : public Curve<Ddata> {
public:
		Dispersion( const std::string& fname )
			: Curve(fname) {}

		float Om( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), Ddata(per) );
			if( itupp == dataV.end() ) {
				if( fabs((*(itupp-1)).per - per) < 0.1 ) return (*(itupp-1)).om;
				else return Ddata::NaN;
			}
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return Ddata::NaN;
			const auto itlow = itupp - 1;
			float odiff = (*itupp).om - (*itlow).om, Tdiff = (*itupp).per - (*itlow).per;
			float deriv = odiff / Tdiff;
			return (*itlow).om + ( per-(*itlow).per ) * deriv;
		}

		float Sdens( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), Ddata(per) );
			if( itupp == dataV.end() ) {
				if( fabs((*(itupp-1)).per - per) < 0.1 ) return (*(itupp-1)).sdensity;
				else return Ddata::NaN;
			}
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return Ddata::NaN;
			const auto itlow = itupp - 1;
			float sdiff = (*itupp).sdensity - (*itlow).sdensity, Tdiff = (*itupp).per - (*itlow).per;
			float deriv = sdiff / Tdiff;
			return (*itlow).sdensity + ( per-(*itlow).per ) * deriv;
		}

		float Deriv_k2om( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), Ddata(per) );
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return Ddata::NaN;
			const auto itlow = itupp - 1;
			return ( (*itupp).k_ang - (*itlow).k_ang ) / ( (*itupp).om - (*itlow).om );
		}

		bool Deriv_k2om( KDeriv& curveout ) {
			curveout.clear();
			if( dataV.size() <= 1 ) return false;
			for( auto it=dataV.begin()+1; it<dataV.end(); it++ ) {
				float per = 0.5 * ( (*it).per + (*(it-1)).per );
				float om = Om( per );
				float deriv = Deriv_k2om( per );
				curveout.push_back(per, deriv, om);
			}
		}

		void ComputeSpectrum( const std::string& sacname ) {
			SacRec sac(sacname);
			sac.Load();
			SacRec sac_am;
			sac.ToAm(sac_am);
			sac_am.Smooth(0.0002, sac);
			float *sigsac = sac.sig.get();
			for( auto& data : dataV ) {
				float freq = 1./data.per;
				int ifreq = (int)floor(freq/sac.shd.delta+0.5);
				data.sdensity = sigsac[ifreq];
				//std::cout<<freq<<" "<<data.sdensity<<std::endl;
			}
		}

		void Write( const std::string& fname ) {
			std::ofstream fout( fname );
			if( ! fout ) {
				std::cerr<<"Error(Dispersion::Write): cannot write to file "<<fname<<std::endl;;
				exit(-2);
			}
			for( const auto& spP : dataV )
				fout<<spP<<"\n";
		}
};

inline bool fExist( const std::string& fname ) { return ( access( fname.c_str(), F_OK ) != -1 ); }

int main(int argc, char* argv[]) {
	if( argc!=5 && argc!=6 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [phv infile (per phv spectral_density[optional])] [grv outfile] [sacout] [distance || sac_in(for spectrum)] [distance (only when sac_in is given for argv[4])]"<<std::endl;
		exit(-1);
	}
	
	// check the type of argument[4]
	bool spin = false; float distin;
	if( fExist( argv[4] ) ) {
		spin = true;
		if( argc == 6 ) distin = atof( argv[5] );
	} else {
		distin = atof( argv[4] );
	}
	// read phase dispersion
	Dispersion disp(argv[1]);
	if( spin )
		disp.ComputeSpectrum(argv[4]);

	// take derivative of wavenumber wrt omiga
	KDeriv kderivs;
	disp.Deriv_k2om( kderivs );
	std::vector<Point> grvV;
	kderivs.Reciprocal( grvV );

	//std::cout<<" per range: "<<disp.xmin()<<" - "<<disp.xmax()<<std::endl;
	std::cout<<" per range: "<<kderivs.xmin()<<" - "<<kderivs.xmax()<<std::endl;

	// synthetic at distance 500km
	float pi = Ddata::pi, oopi = 1./pi;
	float ph_init = pi/4.;

	// sac header
	SacRec sacout;
	auto& shd = sacout.shd;
	if( spin ) {
		sacout.LoadHD(argv[4]);
		if( argc == 6 ) {
			shd.dist = distin;
			shd.stla=-12345.;
			shd.stlo=-12345.;
		}
	} else {
		shd.npts = 3000;
		shd.delta = 1.;
		shd.dist = distin;
		shd.b = 0.;
		shd.e = shd.npts * shd.delta;
	}
	float& delta = shd.delta;
	float& dist = shd.dist;

	// initialize sac
	sacout.sig.reset( new float[shd.npts]() );
	float* sigsac = sacout.sig.get();
	
	// complete integration
	float df = 0.0001, dom = 2*pi*df;
	for( float om=2*pi/disp.xmax(); om<2*pi/disp.xmin(); om+=dom ) {
	//for( float om=2*pi/20.; om<2*pi/17.; om+=dom ) {
		float per = 2.*pi/om, c = disp.Val(per);
		float sdens = disp.Sdens(per);
		//std::cout<<om/(2.*pi)<<" "<<sdens<<std::endl;
		float wavenum = 2.*pi / (per * c);
		for (int it=0; it<shd.npts; it++) {
			float t = it * delta;
			sigsac[it] += sdens * cos(om*t - wavenum*dist + ph_init);
		}
	}
	for (int it=0; it<shd.npts; it++) sigsac[it] *= dom;

	/* stationary phase
	// affected time window
	float tmax = 0., tmin = shd.npts*delta;
	for( Point& grvP : grvV ) {
		float grv = grvP.vel, tcur = dist / grv;
		if( tcur > tmax ) tmax = tcur;
		else if( tcur < tmin ) tmin = tcur;
	}
	int itmin = std::max( 0., (tmin*0.95-50.)/delta );
	int itmax = std::min( (double)(shd.npts), (tmax*1.05+50.)/delta+1 );
	// main loop on time
	for(int it=itmin; it<itmax; it++) {
		float t = it * delta;
		std::vector<float> persV;
		// Search for points where kderiv=t/dist and store periods of these points into persV
		kderivs.FindValue(t/dist, persV);
		if( persV.size() < 1 ) continue;
		//for( const auto& pers : persV ) {
		for( int iper=0; iper<std::min((size_t)3, persV.size()); iper++ ) {
			float& pers = persV[iper];
			float oms = 2.*pi/pers, c = disp.Val(pers);
			float wavenum = 2.*pi / (pers * c);
			float k2deriv_abs = fabs( kderivs.Deriv_2om(pers) );
			std::cerr<<t<<" "<<pers<<std::endl;
			sigsac[it] += sqrt( 2.*oopi / (dist*k2deriv_abs) ) * cos(-oms*t + wavenum*dist + 0.25*pi);
		}
	}
	*/ 
	/* 1st order Taylor approximation on a single period
		float df = 0.2, hdom = pi * df;
	//float perfactor = exp(0.005);	// log(per) + 0.005
	//for( float per=disp.xmin(); per<disp.xmax(); per=1./(1./per-df) ) {
	//for( float per=kderivs.xmin(); per<kderivs.xmax(); per*=perfactor ) {
	//for( float per=4.9; per<5.1; per*=perfactor ) {
		float per = 5.;
		float oms = 2.*pi/per, c = disp.Val(per);
		float wavenum = 2.*pi / (per * c);
		//float kderiv = 1./c + per/(c*c)*disp.Deriv(per);
		//float kderiv = disp.Deriv_k2om(per);
		float kderiv = kderivs.Val(per);
		float k2deriv_abs = fabs( kderivs.Deriv_2om(per) );
		//std::cout<<per<<" "<<kderiv<<" "<<disp.Om(per)<<" "<<k2deriv<<std::endl;
		//grvV.push_back( Ddata(per, 1/kderiv) );
		float ftwindow = 5., tmin = dist * kderiv - per*ftwindow, tmax = dist * kderiv + per*ftwindow;
		int itmin = std::max(0.f, tmin/delta), itmax = std::min((float)(shd.npts), tmax/delta+1);
		for(int it=itmin; it<itmax; it++) {
			float t = it * delta;
			//float Y = 0.5*hdom * (t - kderiv*dist);
			//sigsac[it] += hdom*oopi * sin(Y)/Y * cos(oms*t - wavenum*dist);
			sigsac[it] += sqrt( 2.*oopi / (dist*k2deriv_abs) ) * cos(-oms*t + wavenum*dist + 0.25*pi);
		}
	//}
	*/
	// output group dispersion
	std::ofstream fout(argv[2]);
	if( ! fout ) {
		std::cerr<<"Error(main): cannot write to file "<<argv[2]<<std::endl;;
		exit(-2);
	}
	for( const auto& grvP : grvV )
		fout<<grvP<<std::endl;

	// output spectrum
	if( spin ) {
		std::string fspname = std::string(argv[4]) + "_sp.txt";
		disp.Write( fspname );
	}

	// output sac
	sacout.Write( argv[3] );

	return 0;
}
