#include "SacRec.h"
#include "StaSacs.h"
#include <iostream>

int main (int argc, char *argv[]) {
   if( argc != 8) {
      std::cout<<"Usage: "<<argv[0]<<" [SAC_Z] [SAC_H1] [SAC_H2] [SAC_D] [sac_type (0=dis, 1=vel, 2=acc)] [water-depth(m)] [SAC_Z_out]"<<std::endl;
      exit(-1);
   }
   
	/*
	// check input freqs
   double f1 = atof(argv[3]);
   double f2 = atof(argv[4]);
   double f3 = atof(argv[5]);
   double f4 = atof(argv[6]);
   if( f1<0 || f1>=f2 || f2>=f3 || f3>=f4 ) {
      std::cout<<"Incorrect corner frequencies"<<std::endl;
      exit(0);
   }
	*/

	float Eperl = 11., Eperu = 20.;	// tilt noise is strong in 20 - 50 sec
	StaSacs ss(argv[1], argv[2], argv[3], argv[4], atoi(argv[5]), atof(argv[6]), 0., Eperl, Eperu);
	//float Eperl = 10., Eperu = 12.;
	//auto res = ss.RemoveTiltCompliance("", Eperl, Eperu, 2000.);
	auto res = ss.RemoveTiltCompliance(std::string(argv[7])+"_noise_coherences", 2000.);
	std::cout<<"direction & coh_t0 & coh_c0 & coh_h & coh_1 = "<<res<<"   fcutoff = "<<ss.fcutoffCompliance()<<std::endl;

	ss.Write(argv[7],std::string(argv[7])+"_interm");

   return 0;
}
