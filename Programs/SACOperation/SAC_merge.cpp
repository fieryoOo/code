#include "SacRec.h"
#include <iostream>

int main( int argc, char* argv[] ) {
	if( argc!=4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sac1] [sac2] [sacout]"<<std::endl;
		exit(-1);
	}

	try {
		SacRec sac1(argv[1]), sac2(argv[2]);
		sac1.Load(); sac2.Load();
		sac1.Merge(sac2);
		sac1.Write(argv[3]);
	} catch(...) {
		std::cerr<<"Error(main): exception detected!"<<std::endl;
		return -2;
	}

	return 0;
}
