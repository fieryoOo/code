#include <SacRec.h>
#include <iostream>

   /* method that performs different types of filters:
	 * type = 0: Lowpass cosine -f3~f4_
	 * type = 1: highpass cosine _f1~f2-
	 * type = 2: bandpass cosine _f1~f2-f3~f4_
	 * type = 3: lowpass butterworth -fc=f3~n=f4_
	 * type = 4: highpass butterworth _fc=f1~n=f2-
	 * type = 5: bandpass butterworth _fcL=f1~nL=f2-fcR=f3~nR=f4_
	 * type = 6: gaussian _fc=f2~fhlen=f3_ */
int main( int argc, char* argv[] ) {
   if( argc != 8 ) {
      std::cerr<<"*** --------------------------------------------------- ***"<<"\n";
      std::cerr<<"* type=0: lowpass cosine ( f3 & f4 )"<<"\n";
      std::cerr<<"* type=1: highpass cosine ( f1 & f2 )"<<"\n";
      std::cerr<<"* type=2: bandpass cosine ( f1 - f4 )"<<"\n";
      std::cerr<<"* type=3: lowpass butterworth ( fc=f3 & n=f4 )"<<"\n";
      std::cerr<<"* type=4: highpass butterworth ( fc=f1 & n=f2 )"<<"\n";
      std::cerr<<"* type=5: bandpass butterworth ( fcL=f2 & fcR=f3 & n=f4 )"<<"\n";
      std::cerr<<"* type=6: gaussian ( fc=f2 & fhlen=f3 )"<<"\n";
      std::cerr<<"*** --------------------------------------------------- ***"<<"\n";
      std::cerr<<"Usage: "<<argv[0]<<" [sac_name] [f1] [f2] [f3] [f4] [type(0-6)] [out_sac_name]"<<std::endl;
      exit(-1);
   }

   SacRec sacrec(argv[1]);
   sacrec.Load();

	int type=atoi(argv[6]);
   float f1=atof(argv[2]), f2=atof(argv[3]), f3=atof(argv[4]), f4=atof(argv[5]);
   sacrec.Filter(f1,f2,f3,f4,type);
	//sacrec.BandpassBTWFilt( f2, f3, f4, true );

   sacrec.Write(argv[7]);

   return 0;
}
