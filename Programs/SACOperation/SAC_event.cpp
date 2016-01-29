#include "SacRec.h"
#include <iostream>

int main( int argc, char* argv[] ) {
	if( argc != 7 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sacname] [event (14 digits)] [evlon] [evlat] [tlen (sec)] [outname]"<<std::endl;
		exit(-1);
	}

	try {
		SacRec sac( argv[1] );
		sac.Load();
		std::string event = argv[2];
		sac.ZoomToEvent( argv[2], atof(argv[3]), atof(argv[4]), 0., atof(argv[5]), argv[2] );
		sac.Write( argv[6] );
	} catch(...) {
		std::cerr<<"Error(main): exception detected!"<<std::endl;
		return -2;
	}

	return 0;
}
