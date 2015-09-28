#include "SacRec.h"
#include <iostream>

int main( int argc, char* argv[] ) {
	if( argc<2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sacname] [(optional) field1 field2 field3...]"<<std::endl;
		exit(-1);
	}

	try {
		SacRec sac( argv[1] );
		sac.LoadHD();
		if( argc == 2 ) { 
			std::cout<<sac.shd;
		} else for(int i=2; i<argc; i++) {
			sac.PrintHD(argv[i], std::cout); std::cout<<" ";
		}
		std::cout<<std::endl;
	} catch( std::exception& e ) {
		std::cerr<<e.what()<<std::endl;
		return -2;
	} catch(...) {
		std::cerr<<"Error(main): unknown exception!"<<std::endl;
		return -3;
	}

	return 0;
}
