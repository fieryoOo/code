#include "SacRec.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

int main( int argc, char* argv[] ) {
   if( argc != 4 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [sac_in] [sig_tb] [sig_te]"<<std::endl;
      exit(-1);
   }

   /* load in sac data */
   SacRec sacin(argv[1]);
   if( ! sacin.Load() ) {
      std::cerr<<"Error: invalid sac file name: "<<argv[1]<<std::endl;
      exit(0);
   }

   /* find peak */
   float wb = atof(argv[2]), we = atof(argv[3]);
   if( wb >= we || wb < sacin.shd.b || we > sacin.shd.e ) {
      std::cerr<<"Error: sig window out of range!"<<std::endl;
      exit(0);
   }
   float ftmp, min, max;
   sacin.MinMax( wb, we, ftmp, min, ftmp, max );
   float sigpeak = std::max(fabs(min), fabs(max));

   /* compute rms noise */
   if( we > 0. ) {
      wb = we + 500.;
      we = std::min( wb + 1000.f, sacin.shd.e );
   } else {
      we = wb - 500.;
      wb = std::max( we - 1000.f, sacin.shd.b );
   }
   if( we - wb < 100. ) {
      std::cerr<<"Error: Noise time series not long enough!"<<std::endl;
      exit(0);
   }
   float rms;
   sacin.RMSAvg(wb, we, rms);

   /* print snr */
   std::cout<< sigpeak / rms;// << std::endl;

   return 0;
}
