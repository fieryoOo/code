#include "SacRec.h"
#include <iostream>

int main( int argc, char* argv[] ) {
	if( argc != 5 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sacname] [tb] [te] [outname]"<<std::endl;
		exit(-1);
	}

	try {
		SacRec sac( argv[1] );
		sac.Load();
		float tb = atof(argv[2]), te = atof(argv[3]);
		sac.cut( tb, te );
		sac.Write( argv[4] );
	} catch(...) {
		std::cerr<<"Error(main): exception detected!"<<std::endl;
		return -2;
	}

	return 0;
}
