#include "SacRec.h"
#include "SeedRec.h"
#include <iostream>

int main( int argc, char* argv[] ) {
	// inputs
	if( argc != 3 ) {
		std::cerr<<"Usage: "<<argv[0]<<"[fseed] [sta name]"<<std::endl;
		exit(-1);
	}
	// extract
	SacRec sac;
	SeedRec seed( argv[1] );
	if( ! seedcur.ExtractSac( argv[2], dinfo.chname, dinfo.sps, dinfo.rec_outname,
				dinfo.resp_outname, gapfrac, sac ) ) {
		return 0;
}
