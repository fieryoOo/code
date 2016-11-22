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
   if( argc != 3) {
		// sac_type: 0=dis, 1=vel, 2=acc
      std::cout<<"Usage: "<<argv[0]<<" [sac-list (SAC_Z SAC_H1 SAC_H2 SAC_D sac_type water_depth(m) aziH1 SAC_Z_out)] [output_corrected_H&D_components? 0=no, 1=yes-with-suffix=iter1]"<<std::endl;
      exit(-1);
   }

	// output flag
	int oflag = atoi(argv[2]);
	std::string suffix = oflag==0?"":("_iter"+std::to_string(oflag));

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
		std::stringstream ss(line);	// stream to all strings to prevent funny things happening (eg., 2012.JUN... -> i=2012, s=JUN)
		std::string sacnameZ, sacnameH1, sacnameH2, sacnameD, ssactype, swaterdep, saziH1, outname;
		if( ! (ss >> sacnameZ >> sacnameH1 >> sacnameH2 >> sacnameD >> ssactype >> swaterdep >> saziH1 >> outname) ) {
			std::cerr<<"Warning(main): format error within input line: "<<line<<std::endl;
			continue;
		}
		int sactype = stoi(ssactype); float waterdep = stof(swaterdep), aziH1 = stof(saziH1);

		std::stringstream report(std::ios_base::app|std::ios_base::in|std::ios_base::out);
		report<<"Producing sac file "<<outname<<"\n";
		bool sacWritten = false;
		try {
			// construct StaSacs object
			StaSacs stasac = StaSacs(sacnameZ, sacnameH1, sacnameH2, sacnameD, sactype, waterdep, aziH1, Eperl, Eperu);
			// remove tilt and compliances
			auto res = stasac.RemoveTiltCompliance(outname+"_noise_coherences"+suffix, 2000.);
			report<<"direction & coh_t0 & coh_c0 & coh_h & coh_1 = "<<res<<"   fcutoff = "<<stasac.fcutoffCompliance()<<"\n";
			if( oflag == 0 ) {
				stasac.Write(outname,outname+"_interm"); sacWritten = true;
			} else {
				stasac.Write(oflag); sacWritten = true;
			}
			// compute Rayleigh wave directionality
			float dazi = 5.0;
			//std::vector<std::pair<float, float>> freqRangeV{{0.03,0.039},{0.039,0.051},{0.051, 0.066},{0.066, 0.086},{0.086, 0.11},{0.12, 0.2},{0.2, 0.4}};
			std::vector<std::pair<float, float>> freqRangeV{{0.0177,0.025},{0.025,0.0354},{0.0354, 0.05},{0.05, 0.0707},{0.0707, 0.1},{0.1, 0.141},{0.141,0.2},{0.2,0.283},{0.283,0.4}};
			// compute Rayleigh wave directionality
			RDirect rdCoh, rdAdm;
			//stasac.RayleighDirectionality(dazi, freqRangeV, 2000., 21000., rdCoh, rdAdm);	// 6 hours
			//stasac.RayleighDirectionality(dazi, freqRangeV, 300., 1800., rdCoh, rdAdm);		// 30 minutes
			stasac.RayleighDirectionality(dazi, freqRangeV, 100., 600., rdCoh, rdAdm);			// 10 minutes
			rdCoh.Write(outname+"_rdCoh"+suffix); rdAdm.Write(outname+"_rdAdm"+suffix);
			report<<outname<<"_rdCoh/rdAdm"<<suffix<<" (Rayleigh wave directionalities)"<<std::endl;
		} catch( const std::exception& e ) {
			report<<"Warning(main): StaSacs operation failed ("<<e.what()<<") and no correction made.\n";
			if( ! sacWritten ) {
				auto copyF = [](const std::string &inname, const std::string &outname) {
					std::ifstream fin(inname, std::ios::binary);
					std::ofstream fout(outname, std::ios::binary);
					fout << fin.rdbuf();
				};
				//SacRec sacZ(sacnameZ); sacZ.Load(); sacZ.Write(outname);
				if( oflag == 0 ) {
					copyF(sacnameZ, outname);
				} else {
					copyF(sacnameZ, sacnameZ+suffix);
					copyF(sacnameH1, sacnameH1+suffix);
					copyF(sacnameH2, sacnameH2+suffix);
					copyF(sacnameD, sacnameD+suffix);
				}
			}
		}
		// report 
		#pragma omp critical
		std::cout<<report.str()<<std::endl;
	}

   return 0;
}
