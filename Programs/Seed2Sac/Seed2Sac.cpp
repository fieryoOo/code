#include "SeedRec.h"
#include <iostream>
#include <string>

int main( int argc, char* argv[] ) {
	if( argc != 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sta name]"<<std::endl;
		exit(-1);
	}
	std::string staname(argv[1]);
	SeedRec seedrec("TEST/TA2012_2012.SEP.1.751516.seed");
	SacRec sac;
	float gapfrac;
	std::string rec_outname("TEST/rec.txt"), resp_outname("TEST/RESP_" + staname);
	if( seedrec.ExtractSac(staname, "LHZ", 1, gapfrac, rec_outname, resp_outname, sac) ) {
		sac.Write( std::string("TEST/" + staname + ".SAC").c_str() );
		// remove resp
		sac.RmRESP( resp_outname.c_str(), 5., 80. );
		sac.ZoomToEvent( "20120901000000", -12345., -12345., 1000., 83000. );
		sac.Write( std::string("TEST/ft_" + staname + ".SAC").c_str() );
	}
	return 0;
}
