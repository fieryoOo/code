#include "SacRec.h"
#include "Curve.h"
#include "BSpline.h"
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

		// output files
		std::string outname( argv[1] );
		std::ofstream fphsm(outname + "_phsm");
		//std::ofstream famst(outname + "_amst");
		if( !fphsm ) {//|| !famst ) {
			std::cerr<<"Error(main): cannot write to file "<<outname<<"_phsm"<<std::endl;// or "<<outname<<"_amst"<<std::endl;
			return -2;
		}

		// compute spectrum for each 100sec window
		float hlen = 20., step = 20.;
		float timmin = std::max( (double)step, (double)dist*0.2 );
		float timmax = std::min( (double)sac.shd.e, (double)dist*3. );
		for( float t = timmin; t<=timmax; t+=step ) {
			// apply taper and do FFT
			SacRec saccur(sac); 
			saccur.gauTaper(t, hlen);
			SacRec sac_am, sac_ph;
			saccur.ToAmPh(sac_am, sac_ph);

			// unwrap phase
			sac_ph.Unwrap();

			float fsmhlen = 0.008;
			// store unwrapped phase info into a Curve
			auto& shd = sac_ph.shd;
			float *sigph = sac_ph.sig.get();
			Curve<PointC> curve_ph; curve_ph.reserve(shd.npts);
			for(int i=0; i<shd.npts; i++) {
				curve_ph.push_back( PointC(i*shd.delta, sigph[i]) );
			}
			sigph = nullptr; sac_ph.clear();

			// compute B-Spline on phase
			Curve<PointC> curve_ph_sm;
			BSpline bspline( curve_ph, 1400, 1 );
			bspline.Evaluate( curve_ph_sm );
			//curve_ph.Write("debug_phunwrapped");
			//curve_ph_sm.Write("debug_phspline");

			// estimate rms misfit at each point
			Curve<PointC> curve_diff = curve_ph_sm - curve_ph;
			Curve<PointC> curve_ph_rms = curve_diff.BinRMS(fsmhlen);
			//curve_ph_rms.Write("debug_rms");
			curve_diff.clear();

			// output
			const float vel = dist / t;
			/*
			for( const auto& pc : curve_ph_rms )
				fphsm << 1./pc.x << " " << vel << " " << pc.y << "\n";
			*/
			auto iter = curve_ph_rms.begin();
			sac_am.Transform2( [&](const float f, float& am) {
				fphsm << 1./f << " " << vel << " " << (iter++)->y << " " << am << "\n";
			} );
			fphsm << "\n"; //famst << "\n"; // add blank line for gnuplot pm3d
		}
	} catch(...) {
		std::cerr<<"Error(main): exception detected!"<<std::endl;
		return -3;
	}

	return 0;
}
