#include "SacRec.h"
#include "Dispersion.h"
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <stdexcept>

/* -------------------- the RNG class-------------------- */
class Rand {
   std::default_random_engine generator1;
   std::uniform_real_distribution<float> d_uniform;
   std::normal_distribution<float> d_normal;
public:
   Rand() /* add a true random number from std::random_device to the time seed to ensure thread safety */
      : generator1( std::chrono::system_clock::now().time_since_epoch().count() + std::random_device{}() )
      , d_uniform(0., 1.)
      , d_normal(0., 1.) {}

   //~Rand() {} // should not be defined!!!

   inline float Uniform() { return d_uniform(generator1); }
   inline float Normal() { return d_normal(generator1); }

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
	KDeriv kderivs, grvs;
	disp.Deriv_k2om( kderivs );
	//std::vector<Point> grvV;
	kderivs.Reciprocal( grvs );

	//std::cout<<" per range: "<<disp.xmin()<<" - "<<disp.xmax()<<std::endl;
	std::cout<<" per range: "<<kderivs.xmin()<<" - "<<kderivs.xmax()<<std::endl;

	// synthetic at distance 500km
	float pi = M_PI, oopi = 1./pi;
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
	Rand randO;
	float noiselevel = 0.0;
	for (int it=0; it<shd.npts; it++) {
		sigsac[it] = sigsac[it]*dom + noiselevel * (2.*randO.Normal()-1.);
	}

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
	grvs.Write( argv[2] );
	//std::ofstream fout( argv[2] );
	//for( const auto& grv : grvV ) fout << grv << "\n";

	// output spectrum
	if( spin ) {
		std::string fspname = std::string(argv[4]) + "_sp.txt";
		disp.Write( fspname );
	}

	// output sac
	sacout.Write( argv[3] );

	return 0;
}
