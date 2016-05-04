#include "SacRec.h"
#include "StaSacs.h"
#include <iostream>

int main (int argc, char *argv[]) {
   if( argc != 7) {
      std::cout<<"Usage: "<<argv[0]<<" [SAC_Z] [SAC_D] [sac_type (0=dis, 1=vel, 2=acc)] [water-depth(m)] [zero-pha?(0=no, 1=yes)] [SAC_Z_out]"<<std::endl;
      exit(-1);
   }
   
	bool zeropha = atoi(argv[5])!=0;
	StaSacs ss(argv[1], "", "", argv[2], atoi(argv[3]), atof(argv[4]));
	float Eperl = 15., Eperu = 40.;	// problem: tilt noise is strong in this band too!!
	auto res = ss.RemoveCompliance(std::string(argv[6])+"_noise_coherences", zeropha, Eperl, Eperu, 2000.);
	std::cout<<"coh_c = "<<res<<std::endl;

	ss.Write(argv[6]);

   return 0;
}
