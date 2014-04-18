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
   FocalInfo<float> finfo(strike, dip, rake, dep);
   std::cout<<"### Input Focal info = "<<finfo<<". ###"<<std::endl;

   /* run rad_pattern_r */
   std::vector< std::vector<AziData> > per_azi_pred;
   RadPattern rp;
   rp.Predict( type, argv[2], argv[3], finfo, perlst, per_azi_pred );

   /* output */
   std::ofstream fout(argv[9]);
   if( ! fout ) {
      std::cerr<<"Cannot write to file "<<argv[9]<<std::endl;
      exit(0);
   }
   for(int iper=0; iper<per_azi_pred.size(); iper++) {
      float per_cur = perlst.at(iper);
      std::vector<AziData>& dataV_cur = per_azi_pred.at(iper);
      for(int idat=0; idat<dataV_cur.size(); idat++) {
	 AziData& data_cur = dataV_cur[idat];
	 fout<<data_cur.azi<<"  "<<data_cur.misG<<" "<<data_cur.misP<<" "<<data_cur.A<<"  "<<per_cur<<std::endl;
      }
      fout<<std::endl<<std::endl;
   }
   fout.close();

   return 0;
}
