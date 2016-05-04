#include "SacRec.h"
#include "StaSacs.h"
#include <iostream>

int main (int argc, char *argv[]) {
   if( argc != 9) {
      std::cout<<"Usage: "<<argv[0]<<" [SAC_Z] [SAC_H1] [SAC_H2] [SAC_D] [sac_type (0=dis, 1=vel, 2=acc)] [water-depth(m)] [azi_H1] [SAC_Z_out]"<<std::endl;
      exit(-1);
   }
   
	// remove tilt and compliance
	StaSacs ss(argv[1], argv[2], argv[3], argv[4], atoi(argv[5]), atof(argv[6]), atof(argv[7]));
	float Eperl = 15., Eperu = 40.;	// problem: tilt noise is strong in this band too!!
	//float Eperl = 10., Eperu = 12.;
	//auto res = ss.RemoveTiltCompliance("", Eperl, Eperu, 2000.);
	std::string oname(argv[8]);
	auto res = ss.RemoveTiltCompliance(oname+"_noise_coherences", Eperl, Eperu, 2000.);
	std::cout<<"direction & coh_t0 & coh_c0 & coh_h & coh_1 = "<<res<<"   fcutoff = "<<ss.fcutoffCompliance()<<std::endl;
	ss.Write(oname,oname+"_interm");

	// compute Rayleigh wave directionality
	float dazi = 2.0;
	std::vector<std::pair<float, float>> freqRangeV{{0.03,0.039},{0.039,0.05},{0.05, 0.065},{0.065, 0.085},{0.085, 0.11},{0.12, 0.2},{0.2, 0.4}};
	auto resVV = ss.RayleighDirectionality(dazi, freqRangeV);
	std::ofstream fout(oname+"_RDirect");
	for( const auto &aziLine : resVV ) {
		for( const auto &val : aziLine ) fout<<val<<" ";
		fout<<"\n";
	}
	std::cout<<oname<<"_RDirect (Rayleigh wave directionalities)"<<std::endl;

   return 0;
}
