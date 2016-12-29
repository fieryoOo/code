#include "SacRec.h"
#include "Dispersion.h"
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <stdexcept>

int main(int argc, char* argv[]) {
	if( argc!=3 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [phv infile (per phv)] [grv outfile]"<<std::endl;
		exit(-1);
	}
	
	// take input
	std::string name_grvout(argv[2]);

	// read phase dispersion
	Dispersion disp(argv[1]);

	// take derivative of wavenumber wrt omiga
	KDeriv kderivs, grvs;
	disp.Deriv_k2om( kderivs );

	// predict group
	//std::vector<PointC> grvV;
	kderivs.Reciprocal( grvs );
	grvs.Write( name_grvout );
	return 0;
}
