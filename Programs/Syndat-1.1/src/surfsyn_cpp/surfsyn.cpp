#include "SynGenerator.h"
#include <iostream>
#include <string>
#include <cstdio>

int main( int argc, char* argv[] ) {
	// input params
	if( argc != 6 && argc != 7 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [fparam] [feigen] [wavetype] [mode# (0=fundamental)] [depth] [fix_vel (optional)]"<<std::endl;
		return -1;
	}
	// mode #
	int mode = atoi(argv[4]) + 1;
	if( mode != atof(argv[4])+1 )
		throw std::runtime_error("invalid mode# input (expecting integer)");
	// depth
	float depth = atof(argv[5]);
	if( depth<0 || depth>100 )
		throw std::runtime_error("invalid depth input (expecting 0.0-100.0)");
	if( argc==7 ) {
		float vel = atof(argv[6]);
		if( vel<0.3 || vel>10. ) 
			throw std::runtime_error(std::string("invalid fix_vel input: ")+argv[6]);
	}

	// construct SynGenerator object
	SynGenerator synG( argv[1], argv[2], argv[3][0], mode, depth );

	synG.TraceAllSta();

	synG.ComputeSyn();

	return 0;
}
