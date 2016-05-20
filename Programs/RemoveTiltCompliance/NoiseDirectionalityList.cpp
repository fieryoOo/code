#include "SacRec.h"
#include "StaSacs.h"
#include "MyOMP.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
//#define DEBUG

int main (int argc, char *argv[]) {
   if( argc != 2) {
		// sac_type: 0=dis, 1=vel, 2=acc
      std::cout<<"Usage: "<<argv[0]<<" [sac-list (SAC_Z SAC_H1 SAC_H2 SAC_D sac_type water_depth(m) aziH1 SAC_Z_out)]"<<std::endl;
      exit(-1);
   }

	// read in sac list   
	std::ifstream flst( argv[1] );
	if( ! flst )
		throw std::runtime_error("Error(main): IO failed on file "+std::string(argv[1]));
	std::vector<std::string> flstV;
	for(std::string line; std::getline(flst, line); ) 
		flstV.push_back(std::move(line));

	// loop through the sac list
	float Eperl = 11., Eperu = 20.;
	//for( auto line : flstV ) {
	#pragma omp parallel for schedule(dynamic, 1)
	for(int i=0; i<flstV.size(); i++) {
		// load record from input flist
		const auto &line = flstV[i];
		std::stringstream ss(line);
		std::string sacnameZ, sacnameH1, sacnameH2, sacnameD, outname;
		int sactype; float waterdep, aziH1;
		if( ! (ss >> sacnameZ >> sacnameH1 >> sacnameH2 >> sacnameD >> sactype >> waterdep >> aziH1 >> outname) ) {
			std::cerr<<"Error(main): format error within input line: "<<line<<std::endl;
			continue;
		}
		// remove tilt and compliances
		StaSacs stasac(sacnameZ, sacnameH1, sacnameH2, sacnameD, sactype, waterdep, aziH1, Eperl, Eperu);
		std::stringstream report(std::ios_base::app|std::ios_base::in|std::ios_base::out);
		report<<"Producing sac file "<<outname<<"\n";
		try {
			auto res = stasac.RemoveTiltCompliance(outname+"_noise_coherences", 2000.);
			report<<"direction & coh_t0 & coh_c0 & coh_h & coh_1 = "<<res<<"   fcutoff = "<<stasac.fcutoffCompliance()<<"\n";
		} catch( const std::exception& e ) {
			report<<"Warning(main): rmTiltCompliance failed ("<<e.what()<<") and no correction made.\n";
		}
		stasac.Write(outname,outname+"_interm");
		// compute Rayleigh wave directionality
		float dazi = 5.0;
		//std::vector<std::pair<float, float>> freqRangeV{{0.03,0.039},{0.039,0.051},{0.051, 0.066},{0.066, 0.086},{0.086, 0.11},{0.12, 0.2},{0.2, 0.4}};
		std::vector<std::pair<float, float>> freqRangeV{{0.0177,0.025},{0.025,0.0354},{0.0354, 0.05},{0.05, 0.0707},{0.0707, 0.1},{0.1, 0.141},{0.141,0.2},{0.2,0.283},{0.283,0.4}};
		auto rdirect = stasac.RayleighDirectionality(dazi, freqRangeV, 300., 1800.);
		rdirect.Write(outname+"_RDirect");
		report<<outname<<"_RDirect (Rayleigh wave directionalities)"<<std::endl;
		// report 
		#pragma omp critical
		std::cout<<report.str()<<std::endl;
	}

   return 0;
}
