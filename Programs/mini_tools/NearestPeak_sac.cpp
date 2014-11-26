#include "SacRec.h"
#include <iostream>
#include <cmath>

int main( int argc, char* argv[] ) {
	if( argc != 3 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [fsac] [Tpeak(sec)]"<<std::endl;
		exit(-1);
	}

	try {
		// load
		SacRec sac(argv[1]);
		sac.Load();
		// locate Tpeak
		float delta = sac.shd.delta;
		float Tpeak = atof(argv[2]);
		int ipeak = (int)floor((Tpeak-sac.shd.b)/delta+0.5);
		if( ipeak<0 || ipeak>sac.shd.npts ) {
			std::cerr<<"Error(main): Tpeak out of sac data range!"<<std::endl;
			exit(-2);
		}
		// nearest peak in backward direction
		int i;
		float* sig = sac.sig.get();
		for( i=ipeak-1; i>0; i-- )
			if( sig[i]>sig[i-1] && sig[i]>sig[i+1] ) break;
		float peak_b = sig[i];
		float Tpeak_b = sac.shd.b + i*delta;
		// nearest peak in the forward direction
		for( i=ipeak; i<sac.shd.npts-1; i++ )
			if( sig[i]>sig[i-1] && sig[i]>sig[i+1] ) break;
		float peak_f = sig[i];
		float Tpeak_f = sac.shd.b + i*delta;
		// choose the closer one
		float peak;
		if( (Tpeak-Tpeak_b)<(Tpeak_f-Tpeak) ) {
			Tpeak = Tpeak_b; peak = peak_b;
		} else {
			Tpeak = Tpeak_f; peak = peak_f;
		}
		// output
		std::cout<<Tpeak<<" "<<peak<<std::endl;
	} catch(std::exception& e) {
		std::cerr<<e.what()<<std::endl;
	} catch(...){
		std::cerr<<"Unknown exception!"<<std::endl;
	}

	return 0;
}
