/* testing the methods of the SacRec object */

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
   if( sacrec.sig ) std::cerr<<sacrec.shd.dist<<" "<<sacrec.shd.o<<" "<<sacrec.sig[0]<<" "<<sacrec.sig[shd.npts-1]<<std::endl;

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
   sacrec2.Write("test_shift.SAC");

   /* search for min and max in a given window */
   float tmin, tmax, min, max;
   if( sacrec2.MinMax(9000., 15000., tmin, min, tmax, max) );
      std::cerr<<"9000 - 15000 sec: min="<<min<<" at "<<tmin<<" sec; max="<<max<<" at "<<tmax<<"sec."<<std::endl;

   /* compute rms average in a given window */
   float rms;
   if( sacrec2.RMSAvg(9000., 15000., rms) )
      std::cerr<<"                  rms="<<rms<<std::endl;

   /* test filters */
   // bandpass out-of-place
   SacRec sacrec3;
   sacrec.Filter(0.04, 0.05, 0.1, 0.12, sacrec3);
   sacrec3.Write("bandpass.SAC");
   // bandpass in-place
   sacrec2 = sacrec;
   sacrec2.Filter(0.04, 0.05, 0.1, 0.12);
   sacrec2.Write("bandpass_inplace.SAC");
   // lowpass
   sacrec.Filter(-1., 0., 0.05, 0.06, sacrec3);
   sacrec3.Write("lowpass.SAC");
   // gaussian
   sacrec.Filter(-1., 0.075, 0.025, -1., sacrec3);
   sacrec3.Write("gaussian.SAC");

   return 0;
}
