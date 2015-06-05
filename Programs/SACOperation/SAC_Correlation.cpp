#include "SacRec.h"
#include <iostream>

int main( int argc, char* argv[] ) {
	if( argc!=3 && argc!=5 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sac1] [sac2] [tmin (optional)] [tmax (optional)]"<<std::endl;
		exit(-1);
	}

	try {
		SacRec sac1(argv[1]), sac2(argv[2]);
		sac1.Load(); sac2.Load();
		// decide multiplier
		if( argc == 3 ) {
			std::cout<<sac1.Correlation(sac2)<<std::endl;
		} else {
			std::cout<<sac1.Correlation(sac2, atof(argv[3]), atof(argv[4]))<<std::endl;
		}
		// normalize and output
	} catch(...) {
		std::cerr<<"Error(main): exception detected!"<<std::endl;
		return -2;
	}

	return 0;
}
