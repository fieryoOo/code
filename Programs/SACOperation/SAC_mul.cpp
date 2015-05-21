#include "SacRec.h"
#include <iostream>

int main( int argc, char* argv[] ) {
	if( argc != 4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sacname] [mul] [outname]"<<std::endl;
		exit(-1);
	}

	try {
		SacRec sac( argv[1] );
		sac.Load();
		float mul = atof(argv[2]);
		sac.Mul( mul );
		sac.Write( argv[3] );
	} catch(...) {
		std::cerr<<"Error(main): exception detected!"<<std::endl;
		return -2;
	}

	return 0;
}
