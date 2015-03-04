#include <SacRec.h>
#include <iostream>

int main( int argc, char* argv[] ) {
   if( argc != 7 ) {
      std::cerr<<"*** --------------------------------------------------- ***"<<"\n";
      std::cerr<<"* lowpass when ( (f1==-1. || f2==-1.) && (f3>0. && f4>0.) )"<<"\n";
      std::cerr<<"* bandpass when ( f1>=0. && f2>0. && f3>0. && f4>0. )"<<"\n";
      std::cerr<<"* gaussian when ( f1==-1. && f4==-1. ) and "<<"\n"
	       <<"  f2 = center freqency and f3 = frequency half length"<<"\n";
      std::cerr<<"*** --------------------------------------------------- ***"<<"\n";
      std::cerr<<"Usage: "<<argv[0]<<" [sac_name] [f1] [f2] [f3] [f4] [out_sac_name]"<<std::endl;
      exit(-1);
   }

   SacRec sacrec(argv[1]);
   sacrec.Load();

   float f1=atof(argv[2]), f2=atof(argv[3]), f3=atof(argv[4]), f4=atof(argv[5]);
   sacrec.Filter(f1,f2,f3,f4);

   sacrec.Write(argv[6]);

   return 0;
}
