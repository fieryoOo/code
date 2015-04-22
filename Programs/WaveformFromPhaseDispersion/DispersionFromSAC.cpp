#include "SacRec.h"
#include "Dispersion.h"
#include "BSpline.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

inline int nint( float f ) { return (int)floor(f+0.5); }
int main(int argc, char* argv[]) {
	if( argc != 4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sac_in] [tmax] [permax]\n";
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
		sigsac = nullptr;
//sac.Write(std::string(argv[1])+"_tapered");

		// FFT
		SacRec sac_am, sac_am_o, sac_ph_o, sac_ph;
		sac.ToAmPh( sac_am_o, sac_ph_o );
sac_am_o.Write( std::string(argv[1]) + ".am" );
sac_ph_o.Write( std::string(argv[1]) + ".ph" );
		sac_ph = sac_ph_o;
		//sac.clear(); sac_am_o.clear(); 
		sac.sig.reset();

		// unwrap phase
		float fsmhlen = 0.008;
		float twopi = 2. * M_PI;
		auto& shd = sac_ph.shd;
		float *sigph = sac_ph.sig.get();
		for(int i=1; i<shd.npts; i++) {
			sigph[i] = sigph[i] - nint( (sigph[i]-sigph[i-1])/twopi ) * twopi;
		}
//sac_ph.Write("Phase.SAC");
		//SacRec sac_tmp;
		//sac_ph.Smooth(fsmhlen, sac_tmp);
		//sac_ph = std::move(sac_tmp);

		// store phase info into a Curve
		Curve<Point> curve_ph; curve_ph.reserve(shd.npts);
		for(int i=0; i<shd.npts; i++) {
			curve_ph.push_back(i*shd.delta, sigph[i]);
		}
		sigph = nullptr; sac_ph.clear();
		// compute B-Spline on phase
		Curve<Point> curve_ph_sm;
		BSpline bspline( curve_ph, 1200, 1 );
		bspline.Evaluate( curve_ph_sm );

curve_ph.Write("Phase.SAC_unwrapped");
curve_ph_sm.Write("Phase.SAC_spline");
		// estimate rms misfit at each point
		curve_ph_sm -= curve_ph;
		//curve_ph_sm = curve_ph - curve_ph_sm;
		Curve<Point> curve_ph_rms;
		curve_ph_sm.BinRMS( fsmhlen, curve_ph_rms );
curve_ph_rms.Write("Phase.SAC_rms");
		curve_ph_sm.clear();

		// convert phase to velocity, discarding points with misfit > mismax
		float dist = sac.shd.dist;
		float mismax = 0.2;								// maximum allowed rms misfit in phase
		float oo2pi = 1./twopi, vel = 3.8;			// reference velocity
		auto iterbeg = curve_ph.lower_bound(1./atof(argv[3]));							// ignore per > 100.
		float n2pi = nint( dist*(iterbeg->x)/vel - ((iterbeg->y)*oo2pi+0.125) );	// pick the branch closet to ref vel
		Dispersion disp;
		for(auto iter=iterbeg; iter<curve_ph.end(); iter++) {
			auto iterrms = curve_ph_rms.begin() + (iter-curve_ph.begin());
			if( iterrms->y > mismax ) continue;
			float &freq = iter->x;
			float ftmp1 = dist*freq, ftmp2 = (iter->y)*oo2pi + 0.125;	// pi/4 shift for AN
			vel = ftmp1 / (ftmp2 + n2pi);
			disp.push_back(1./freq, vel);
		}
		curve_ph.clear();
		curve_ph_rms.clear();

		// smooth amplitude ( for invalidating noise in the next step )
		//sac_am_o.Smooth( fsmhlen, sac_am );
		//float *sigam = sac_am.sig.get(), *sigamo = sac_am_o.sig.get();
		
/*
		// work on the phase
		Dispersion disp;
		sigph = sac_ph.sig.get();
		float oo2pi = 1./twopi, vel = 4.5;			// reference velocity
		int ibeg = nint(0.01/shd.delta);				// ignore per > 100.
		float n2pi = nint( dist*ibeg*shd.delta/vel - (sigph[ibeg]*oo2pi+0.125) );	// pick the branch closet to ref vel
		for(int i=ibeg; i<shd.npts; i++) {
			//if( fabs(sigamo[i]-sigam[i]) > sigam[i]*0.05 ) continue;
			float freq = i*shd.delta;
			float ftmp1 = dist*freq, ftmp2 = sigph[i]*oo2pi + 0.125;	// pi/4 shift for AN
			//float n2pi = nint( ftmp1/vel - ftmp2 );
			vel = ftmp1 / (ftmp2 + n2pi);
			disp.push_back(1./freq, vel);
		}
		sigph = nullptr; //sigam = sigamo = nullptr;
*/
		disp.Sort();
/*
		BSpline bspline2( disp, 1200, 1 );
		bspline2.Evaluate( curve_ph_sm );
curve_ph_sm.Write("aa");
*/


		// take derivative of wavenumber wrt omiga
		KDeriv kderivs, grvs;
		disp.Deriv_k2om( kderivs );
		kderivs.Reciprocal( grvs );

		// check group curve and remove noisy points
/*
		float maxVerror = 0.01;
		Point grvP;
		for( grvs.rewind(); grvs.get(grvP); grvs.next() ) {
			float &per = grvP.x, &grv = grvP.y;
			if( per < 1.5/sac_am_o.shd.e ) continue;
//if( grv<0.5 || grv>5. ) continue;
			SacRec sac_am = sac_am_o, sac_env;
			sac_am.gauTaper( 1./per, fsmhlen );		// apply gaussian filter
			sac_env.shd = sac.shd;						// use header of the original sac file
			sac_env.FromAmPh( sac_am, sac_ph_o, 2 );	// and get envelope
			float tgr = dist/grv;	// allow 5% error in grv
			if( tgr<10. || tgr>sac_env.shd.e-10 ) continue;
			float tmin, tmax, ftmp;
			sac_env.MinMax(tgr*(1-2.*maxVerror), tgr*(1+2.*maxVerror), tmin, ftmp, tmax, ftmp);	// search +/- maxVerror percent for local maximum
			//if( tmax<tgr*(1-maxVerror) || tmax>tgr*(1+maxVerror) ) continue;	// local maximum not found in range

			std::cout<<per<<" "<<grv<<" "<<disp.Val(per)<<"\n";
		}
*/
		// output
		std::string outname( argv[1] );
		disp.Write( outname + "_phv" );
		grvs.Write( outname + "_grv" );

	} catch( std::exception& e ) {
		std::cerr<<"Error(main): "<<e.what()<<"\n";
	} catch(...) {
		std::cerr<<"Error(main): unknown exception detected!\n";
	}

	return 0;
}
