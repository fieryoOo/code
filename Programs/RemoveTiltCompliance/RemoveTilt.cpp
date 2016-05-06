#include "SacRec.h"
#include "StaSacs.h"
#include <iostream>

int main (int argc, char *argv[]) {
   if( argc != 7) {
      std::cout<<"Usage: "<<argv[0]<<" [SAC_Z] [SAC_H1] [SAC_H2] [sac_type (0=dis, 1=vel, 2=acc)] [zero-phase?(0=no, 1=yes)] [SAC_Z_out]"<<std::endl;
      exit(-1);
   }
   
	bool zeropha = atoi(argv[5])!=0;
	float Eperl = 11., Eperu = 20.;	
	StaSacs ss(argv[1], argv[2], argv[3], "", atoi(argv[4]), 0, 0., Eperl, Eperu);
	auto res = ss.RemoveTilt(std::string(argv[6])+"_noise_coherences", 2000.);
	std::cout<<"direction & coh_t = "<<res.x<<" "<<res.y<<std::endl;

	ss.Write(argv[6]);

   return 0;
}
