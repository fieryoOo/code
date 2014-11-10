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

   /* reshape and shift the signal into sacrec2 */
   SacRec sacrec2( sacrec1 );
   double pi = 4.*atan(1.), mper = pi*sacrec2.shd.delta/10000.;
   for(int i=0; i<sacrec2.shd.npts; i++) sacrec2.sig[i] *= sin(i*mper);
   sacrec2.shd.nzmsec += 100000*1000; // shift by 5000 sec
   sacrec2.UpdateTime();

   /* merge and write */
   //sacrec1.Merge(sacrec2);
   //sacrec1.Write("merged1.SAC");
   sacrec1.merge(sacrec2);
   int Nholes = sacrec1.arrange("merged_rec1");
   std::cerr<<Nholes<<std::endl;
   sacrec1.Write("merged2.SAC");

   return 0;
}
