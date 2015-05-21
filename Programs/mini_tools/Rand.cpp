#include "Rand.h"
#include <iostream>

int main( int argc, char* argv[] ) {
	if( argc != 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [0 (uniform) or 1 (normal)]\n";
		return -1;
	}

	int rtype=atoi(argv[1]);

	Rand randO;
	switch( rtype ) {
		case 0:
			std::cout<<randO.Uniform()<<"\n";
			break;
		case 1:
			std::cout<<randO.Normal()<<"\n";
			break;
		default:
			std::cerr<<"Error(main): unknown type rtype="<<rtype<<"\n";
			return -2;
	}

	return 0;
}
