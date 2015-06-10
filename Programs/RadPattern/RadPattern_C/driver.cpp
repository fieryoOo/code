#include "RadPattern.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>


inline int nint(float datain) { return (int)floor(datain+0.5); }

int main( int argc, char* argv[] ) {
   if( argc != 10 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [R/L] [eigen_file (.R/L)] [phvel_file (.R/L.phv)]";
      std::cerr<<" [per_lst] [strike] [dip] [rake] [depth] [out_name]"<<std::endl;
      exit(-1);
   }

   /* read in type */
   char type = argv[1][0];
   if( type != 'R' && type != 'L' ) {
      std::cerr<<"Unknown type: "<<type<<std::endl;
      exit(0);
   }

   /* read in per.lst */
   std::vector<float> perlst;
   std::ifstream fin(argv[4]);
   if( ! fin ) {
      std::cerr<<"Error(main): Cannot read from file "<<argv[4]<<std::endl;
      exit(0);
   }
   for(std::string line; std::getline(fin, line); ) {
      float pertmp;
      sscanf(line.c_str(), "%f", &pertmp);
      perlst.push_back(pertmp);
   }
   fin.close();
   std::cout<<"### "<<perlst.size()<<" periods read in. ###"<<std::endl;

   /* read in focal info */
   float strike = atof(argv[5]), dip = atof(argv[6]), rake = atof(argv[7]), dep = atof(argv[8]);
   /*
   int strike = nint(strikein), dip = nint(dipin), rake = nint(rakein), dep = nint(depin);
   if( strike!=strikein || dip!=dipin || rake!=rakein || dep!=depin ) {
      std::cerr<<"Warning(main): integer expected for strike/dip/rake/depth. Corrected to the nearest integer(s)!"<<std::endl;
   }
   */
   //FocalInfo<float> finfo(strike, dip, rake, dep);
   std::cout<<"### Input Focal info = ("<<strike<<" "<<dip<<" "<<rake<<" "<<dep<<"). ###"<<std::endl;

	//bool GetPred( const float per, const float azi,	float& grt, float& pht, float& amp ) const;
   /* run rad_pattern_r */
   RadPattern rp;
   rp.Predict( type, argv[2], argv[3], strike, dip, rake, dep, perlst );
	rp.OutputPreds( argv[9] );

   return 0;
}
