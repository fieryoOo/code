#include "SacRec.h"
#include <iostream>

int main( int argc, char* argv[] ) {
	if( argc != 3 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sacname] [outname]"<<std::endl;
		exit(-1);
	}

	try {
		SacRec sac( argv[1] );
		sac.Load();
		sac.Envelope();
		sac.Write( argv[2] );
	} catch(...) {
		std::cerr<<"Error(main): exception detected!"<<std::endl;
		return -2;
	}

	return 0;
}
