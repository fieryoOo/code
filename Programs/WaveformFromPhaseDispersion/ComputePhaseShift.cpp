#include "SacRec.h"
#include "Curve.h"
#include <iostream>
#include <fstream>


int main( int argc, char* argv[] ) {
	if( argc != 4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sac1] [sac2] [fout]"<<std::endl;
		return -1;
	}

	// correlate
	SacRec sac1, sac2, saccor;
	sac1.Load(argv[1]);
	sac2.Load(argv[2]);
	sac1.Correlate(sac2, saccor);

	// file for output 
	std::ofstream fout(argv[3]);
	if( ! fout ) {
		std::cerr<<"Error(main): cannot write to file "<<argv[3]<<std::endl;
		return -2;
	}

	// find T of maximum correlation at each period
	float perfactor = 1.05;
	float shdb = saccor.shd.b, delta = saccor.shd.delta;
	for( float per=5.; per<50.; per*=perfactor ) {
		// gauss filt
		const float freq = 1. / per;
		const float dfreq = 1./exp(log(per)-0.15) - freq;
		SacRec sacflt;
		saccor.GaussianFilt( freq, dfreq, sacflt );
		// find maximum
		int imin, imax;
		sacflt.MinMax(imin, imax, -per*0.51, per*0.51);
		if( imax<=0 || imax> sacflt.shd.npts-2 ) continue;
		// the max, left, and right points
		float *sigsac = sacflt.sig.get();
		Point p1( shdb + (imax-1)*delta, sigsac[imax-1] );
		Point p2( shdb + imax*delta, sigsac[imax] );
		Point p3( shdb + (imax+1)*delta, sigsac[imax+1] );
		Parabola pr(p1, p2, p3);
		fout<<per<<" "<<pr.Vertex()<<std::endl;
	}
	
}
