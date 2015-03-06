#include "SacRec.h"
#include "Curve.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

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
		sigsac = nullptr;

		// FFT
		SacRec sac_am, sac_am_o, sac_ph_o, sac_ph;
		sac.ToAmPh( sac_am_o, sac_ph_o );
		sac_ph = sac_ph_o;
		float dist = sac.shd.dist;
		//sac.clear(); sac_am_o.clear(); 
		sac.sig.reset();

		// correct 2pi in phase
		float fsmhlen = 0.008;
		float twopi = 2. * M_PI;
		auto& shd = sac_ph.shd;
		float *sigph = sac_ph.sig.get();
		for(int i=1; i<shd.npts; i++) {
			sigph[i] = sigph[i] - nint( (sigph[i]-sigph[i-1])/twopi ) * twopi;
		}
		sigph = nullptr;
		SacRec sac_tmp;
		sac_ph.Smooth(fsmhlen, sac_tmp);
		sac_ph = std::move(sac_tmp);

		// smooth amplitude ( for invalidating noise in the next step )
		sac_am_o.Smooth( fsmhlen, sac_am );
		float *sigam = sac_am.sig.get(), *sigamo = sac_am_o.sig.get();
		
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
		sigph = sigam = sigamo = nullptr;
		disp.Sort();

		// take derivative of wavenumber wrt omiga
		KDeriv kderivs, grvs;
		disp.Deriv_k2om( kderivs );
		kderivs.Reciprocal( grvs );

		// check group curve and remove noisy points
		float maxVerror = 0.01;
		Point grvP;
		for( grvs.rewind(); grvs.get(grvP); grvs.next() ) {
			float &per = grvP.per, &grv = grvP.vel;
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
