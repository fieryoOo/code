#include "SacRec.h"
#include <iostream>
#include <fstream>
#include <cmath>

inline int nint( float f ) { return (int)floor(f+0.5); }
int main( int argc, char* argv[] ) {
	if( argc != 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sacname]"<<std::endl;
		exit(-1);
	}

	try {
		// load sac
		SacRec sac( argv[1] );
		sac.Load();
		float dist = sac.shd.dist;
		if( dist == SacRec::NaN ) {
			std::cerr<<"Error(main): no distance info in SAC header"<<std::endl;
			return -2;
		}

		// FFT
		SacRec sac_am, sac_ph;
		sac.ToAmPh(sac_am, sac_ph);
		sac.sig.reset(); //sac.clear();

		// output files
		std::string outname( argv[1] );
		std::ofstream fenv(outname + "_env");
		std::ofstream fpha(outname + "_pha");
		if( !fenv || !fpha ) {
			std::cerr<<"Error(main): cannot write to file "<<outname<<"_env or "<<outname<<"_pha"<<std::endl;
			return -2;
		}

		// apply gaussian filter at each period and IFFT to get FTAN diagram
		float twopi = M_PI * 2.;
		for( float per=3; per<100.; per*=1.05 ) {
		//for( float per=3; per<50.; per+=0.5 ) {
         SacRec sac_am_cur = sac_am;
         sac_am_cur.gauTaper( 1./per, 0.008 );		// apply gaussian filter
			SacRec sac_env, sac_pha;
         sac_env.shd = sac_pha.shd = sac.shd;		// use header of the original sac file
         sac_env.FromAmPh( sac_am_cur, sac_ph, 2 );	// get envelope
			sac_pha.FromAmPh( sac_am_cur, sac_ph, 3 );	// and phase (optimization needed!)
			// correct 2 pi in phase
			float *sigpha = sac_pha.sig.get();
			//for(int i=1; i<sac_pha.shd.npts; i++) 
			//	sigpha[i] = sigpha[i] - nint( (sigpha[i]-sigpha[i-1])/twopi ) * twopi;
			// output FTAN amp and phase
			float disodt = dist / sac_env.shd.delta;
			float *sigenv = sac_env.sig.get();
			float velold = 20.1;
			// start at vel=20. and loop down to vel=0.
			for( int i=(int)ceil(disodt/(velold-0.1)); i<sac_env.shd.npts; i++ ) {
				float vel = disodt / i;
				if( velold-vel<0.01 ) continue;	// skip if delta(v) is too small
				//float wavenum = i / (disodt*per);
				fenv << per << " " << vel << " " << sigenv[i] << "\n";
				fpha << per << " " << vel << " " << sigpha[i] << "\n";
				velold = vel;
			}
			fenv << "\n"; fpha << "\n"; // add blank line for gnuplot pm3d
		}
	} catch(...) {
		std::cerr<<"Error(main): exception detected!"<<std::endl;
		return -3;
	}

	return 0;
}
