#include "SacRec.h"
#include "StaSacs.h"
#include <iostream>
//#define DEBUG

int main (int argc, char *argv[]) {
   if( argc != 6) {
      std::cout<<"Usage: "<<argv[0]<<" [SAC_Z] [SAC_H1] [SAC_H2] [sac_type (0=dis, 1=vel, 2=acc)] [SAC_Z_out]"<<std::endl;
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

	StaSacs ss(argv[1], argv[2], argv[3], "", atoi(argv[4]));
	float Eperl = 15., Eperu = 40.;	// problem: tilt noise is strong in this band too!!
	//float Eperl = 10., Eperu = 12.;
	ss.RemoveTilts(Eperl, Eperu, 2000.);

	ss.Write(argv[5]);

   return 0;
}
