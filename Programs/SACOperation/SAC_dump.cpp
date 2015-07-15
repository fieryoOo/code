#include "SacRec.h"
#include <iostream>

int main( int argc, char* argv[] ) {
	if( argc!=2 && argc!=3 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sacname] [outname (optional)]"<<std::endl;
		exit(-1);
	}

	try {
		SacRec sac( argv[1] );
		sac.Load();
		switch( argc ) {
			case 2: sac.Dump(); break;
			case 3: sac.Dump( argv[2] );
		}
	} catch( std::exception& e ) {
		std::cerr<<e.what()<<std::endl;
		return -2;
	} catch(...) {
		std::cerr<<"Error(main): unknown exception!"<<std::endl;
		return -3;
	}

	return 0;
}
