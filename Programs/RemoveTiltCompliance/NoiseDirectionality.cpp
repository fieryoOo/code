#include "SacRec.h"
#include "StaSacs.h"
#include <iostream>

int main (int argc, char *argv[]) {
   if( argc != 9) {
      std::cout<<"Usage: "<<argv[0]<<" [SAC_Z] [SAC_H1] [SAC_H2] [SAC_D] [sac_type (0=dis, 1=vel, 2=acc)] [water-depth(m)] [azi_H1] [SAC_Z_out]"<<std::endl;
      exit(-1);
   }
   
	// remove tilt and compliance
	float Eperl = 11., Eperu = 20.;	// tilt noise is strong in between 20 and 50 sec
	StaSacs ss(argv[1], argv[2], argv[3], argv[4], atoi(argv[5]), atof(argv[6]), atof(argv[7]), Eperl, Eperu);
	std::string oname(argv[8]);
/*
	auto res = ss.RemoveTiltCompliance(oname+"_noise_coherences", 2000.);
	std::cout<<"direction & coh_t0 & coh_c0 & coh_h & coh_1 = "<<res<<"   fcutoff = "<<ss.fcutoffCompliance()<<std::endl;
	ss.Write(oname,oname+"_interm");
	//ss.Write(oname,oname+"_interm",std::string(argv[2])+"_c",std::string(argv[3])+"_c",std::string(argv[4])+"_c");
*/

	// compute Rayleigh wave directionality
	float dazi = 5.0;
	//std::vector<std::pair<float, float>> freqRangeV{{0.03,0.039},{0.039,0.05},{0.05, 0.065},{0.065, 0.085},{0.085, 0.11},{0.12, 0.2},{0.2, 0.4}};
	std::vector<std::pair<float, float>> freqRangeV{{0.0177,0.025},{0.025,0.0354},{0.0354, 0.05},{0.05, 0.0707},{0.0707, 0.1},{0.1, 0.141},{0.141,0.2},{0.2,0.283},{0.283,0.4}};
	RDirect rdCoh, rdAdm;
	ss.RayleighDirectionality(dazi, freqRangeV, 300., 1800., rdCoh, rdAdm);
	//auto rdirect = ss.RayleighDirectionality(dazi, freqRangeV, 450., 3600.);
	rdCoh.Write(oname+"_rdCoh"); rdAdm.Write(oname+"_rdAdm");
	/*
	std::ofstream fout(oname+"_RDirect");
	for( const auto &aziLine : resVV ) {
		for( const auto &val : aziLine ) fout<<val<<" ";
		fout<<"\n";
	}
	*/
	std::cout<<oname<<"_RDirect (Rayleigh wave directionalities)"<<std::endl;

   return 0;
}
