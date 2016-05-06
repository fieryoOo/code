#include "SacRec.h"
#include "StaSacs.h"
#include <iostream>

int main (int argc, char *argv[]) {
   if( argc!=7 && argc!=9 ) {
      std::cout<<"Usage: "<<argv[0]<<" [SAC_Z] [SAC_H1] [SAC_H2] [SAC_D] [sac_type (0=dis, 1=vel, 2=acc)] [water-depth(m)] [tb(optional)] [te(optional)]"<<std::endl;
      exit(-1);
   }
   
	StaSacs ss(argv[1], argv[2], argv[3], argv[4], atoi(argv[5]), atof(argv[6]));
	float twin = 300.;
	if( argc == 7 ) ss.test(twin, 25);
	else ss.test(twin, 25, atof(argv[7]), atof(argv[8]));

   return 0;
}
