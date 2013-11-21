#include "SacRec.h"
#include <cmath>

int main( int argc, char* argv[] )
{
   /* initialize */
   SacRec sacrec(argv[1]);

   /* load header */
   sacrec.LoadHD();
   SAC_HD& shd = sacrec.shd;
   std::cerr<<shd.delta<<" "<<shd.npts<<std::endl;

   /* load signal */
   sacrec.Load();
   if( sacrec.sig ) std::cerr<<shd.dist<<" "<<shd.o<<" "<<sacrec.sig[0]<<" "<<sacrec.sig[shd.npts-1]<<std::endl;

   /* copy constructor */
   SacRec sacrec2( sacrec );
   if( sacrec2.sig ) std::cerr<<sacrec2.shd.dist<<" "<<sacrec2.shd.o<<" "<<sacrec2.sig[0]<<" "<<sacrec2.sig[shd.npts-1]<<std::endl;

   /* assignment operator */
   sacrec2.sig[0] = -1.; sacrec2.sig[shd.npts-1] = -2.; sacrec2.shd.dist = 999.;
   sacrec2 = sacrec;
   if( sacrec2.sig ) std::cerr<<sacrec2.shd.dist<<" "<<sacrec2.shd.o<<" "<<sacrec2.sig[0]<<" "<<sacrec2.sig[shd.npts-1]<<std::endl;

   /* manipulate signal/header and write to file */
   double pi = 4.*atan(1.), mper = pi*shd.delta/10000.;
   for(int i=0; i<shd.npts; i++) sacrec2.sig[i] *= sin(i*mper);
   sacrec2.shd.nzmsec += 5000*1000; // shift by 5000 sec
   sacrec2.Write("temp.SAC");

   return 0;
}
