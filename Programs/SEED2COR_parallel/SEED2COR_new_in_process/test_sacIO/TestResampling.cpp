/* testing the methods of the SacRec object */

#include "SacRec.h"
#include <cmath>

int main( int argc, char* argv[] )
{
   if( argc != 2 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [input_sac]"<<std::endl;
      exit(-1);
   }
   /* read the sac file into sacrec1 */
   SacRec sacrec1(argv[1]);
   sacrec1.Load();

   /* resample to 2 sps */
   sacrec1.Resample(2);

   sacrec1.Write("Resampled_2.SAC");

   return 0;
}
