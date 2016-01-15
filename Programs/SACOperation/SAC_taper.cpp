#include "SacRec.h"
#include <iostream>

int main( int argc, char* argv[] ) {
	if( argc!=7 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sacname] [x1] [x2] [x3] [x4] [outname]"<<std::endl;
		exit(-1);
	}

	float x1 = atof(argv[2]), x2 = atof(argv[3]);
	float x3 = atof(argv[4]), x4 = atof(argv[5]);
	try {
		SacRec sac( argv[1] );
		sac.Load();
		sac.cosTaperL( x1, x2 );
		sac.cosTaperR( x3, x4 );
		sac.Write( argv[6] );
	} catch( std::exception& e ) {
		std::cerr<<e.what()<<std::endl;
		return -2;
	} catch(...) {
		std::cerr<<"Error(main): unknown exception!"<<std::endl;
		return -3;
	}

	return 0;
}
