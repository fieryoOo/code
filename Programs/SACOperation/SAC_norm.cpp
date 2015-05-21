#include "SacRec.h"
#include <iostream>

int main( int argc, char* argv[] ) {
	if( argc!=3 && argc!=5 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sacname] [outname] [tmin (optional)] [tmax (optional)]"<<std::endl;
		exit(-1);
	}

	try {
		SacRec sac( argv[1] );
		sac.Load();
		// decide multiplier
		float mul;
		if( argc == 3 ) {
			mul = 1. / sac.shd.depmax;
		} else {
			float tmin, min, tmax, max;
			sac.MinMax( atof(argv[3]), atof(argv[4]), tmin, min, tmax, max);
			mul = 1. / std::max(max, -min);
			//std::cerr<<" in SAC_norm: "<<tmin<<" "<<min<<"   "<<tmax<<" "<<max<<"\n";
		}
		// normalize and output
		sac.Mul( mul );
		sac.Write( argv[2] );
	} catch(...) {
		std::cerr<<"Error(main): exception detected!"<<std::endl;
		return -2;
	}

	return 0;
}
