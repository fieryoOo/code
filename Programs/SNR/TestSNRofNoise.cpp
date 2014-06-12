#include "SacRec.h"
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>


static constexpr float pi = 3.1415927;

/* --- main --- */
inline int nint( float x ) { return (int)floor(x+0.5); }
int main() {
   std::vector<std::string> saclst = { "TEST/SACs/COR_M12A_P21A.SAC", "TEST/SACs/COR_M12A_S04C.SAC", "TEST/SACs/COR_M12A_W23A.SAC", "TEST/SACs/COR_M12A_Z20A.SAC", "TEST/SACs_OBS/COR_CMB_E04D.SAC", "TEST/SACs_OBS/COR_FN01A_FN12A.SAC", "TEST/SACs_OBS/COR_FN05A_G30A.SAC", "TEST/SACs_OBS/COR_FN07A_J01D.SAC", "TEST/SACs_OBS/COR_FN08A_J30A.SAC", "TEST/SACs_OBS/COR_FN12A_J42A.SAC", "TEST/SACs_OBS/COR_FN14A_J53A.SAC", "TEST/SACs_OBS/COR_FN16A_K02D.SAC", "TEST/SACs_OBS/COR_FN18A_N02D.SAC", "TEST/SACs_OBS/COR_G03D_D04D.SAC", "TEST/SACs_OBS/COR_G30A_N02D.SAC", "TEST/SACs_OBS/COR_I04A_A04D.SAC", "TEST/SACs_OBS/COR_J05D_NLWA.SAC", "TEST/SACs_OBS/COR_J23A_G05D.SAC", "TEST/SACs_OBS/COR_J25A_J34A.SAC", "TEST/SACs_OBS/COR_J26A_J54A.SAC", "TEST/SACs_OBS/COR_J28A_M06A.SAC", "TEST/SACs_OBS/COR_J30A_D04D.SAC", "TEST/SACs_OBS/COR_J31A_J39A.SAC", "TEST/SACs_OBS/COR_J33A_J73A.SAC", "TEST/SACs_OBS/COR_J35A_B05D.SAC", "TEST/SACs_OBS/COR_J36A_J42A.SAC", "TEST/SACs_OBS/COR_J37A_M02A.SAC", "TEST/SACs_OBS/COR_J39A_G03D.SAC", "TEST/SACs_OBS/COR_J41A_J63A.SAC", "TEST/SACs_OBS/COR_J43A_C06D.SAC", "TEST/SACs_OBS/COR_J44A_J61A.SAC", "TEST/SACs_OBS/COR_J46A_D03D.SAC", "TEST/SACs_OBS/COR_J47A_K04D.SAC", "TEST/SACs_OBS/COR_J49A_I02D.SAC", "TEST/SACs_OBS/COR_J50A_PINE.SAC", "TEST/SACs_OBS/COR_J52A_K04D.SAC", "TEST/SACs_OBS/COR_J54A_J04D.SAC", "TEST/SACs_OBS/COR_J57A_F05D.SAC", "TEST/SACs_OBS/COR_J59A_BMO.SAC", "TEST/SACs_OBS/COR_J63A_A04D.SAC", "TEST/SACs_OBS/COR_J67A_B05D.SAC", "TEST/SACs_OBS/COR_J73A_D03D.SAC", "TEST/SACs_OBS/COR_L02D_HAWA.SAC", "TEST/SACs_OBS/COR_M02A_A04D.SAC", "TEST/SACs_OBS/COR_M03A_M02C.SAC", "TEST/SACs_OBS/COR_M07A_HAWA.SAC", "TEST/SACs_OBS/COR_MCCM_PINE.SAC", "TEST/SACs_OBS/COR_PINE_NLWA.SAC" };
   float fhlen = 0.01;
   float fl = 0.03, fh = 0.21, fstep=0.02;
   int nf = nint( (fh - fl) / fstep ) + 1;
   std::ofstream fout("SNR_of_Noise.txt");
   for(int isac=0; isac<saclst.size(); isac++) {
      // load in sac file
      SacRec sacin(saclst[isac].c_str());
      sacin.Load();
      // g filter to different periods
      for(int ifreq=0; ifreq<nf; ifreq++) {
	 float freq = fl + ifreq * fstep;
	 SacRec sacft;
	 sacin.Filter(-1., freq, fhlen, -1., sacft);
	 // compute peak on pos lag
	 float min, max, ftmp;
	 sacft.MinMax(1500., 2500., ftmp, min, ftmp, max);
	 float peak_pos = std::max(fabs(min), fabs(max));
	 // compute peak on neg lag
	 sacft.MinMax(-2500., -1500., ftmp, min, ftmp, max);
	 float peak_neg = std::max(fabs(min), fabs(max));
	 // compute rms on each lag
	 float rms_pos, rms_neg;
	 sacft.RMSAvg(1500., 2500., 1, rms_pos);
	 sacft.RMSAvg(-2500., -1500., 1, rms_neg);
	 fout<<freq<<" "<<peak_pos/rms_pos<<" pos "<<saclst[isac]<<"\n";
	 fout<<freq<<" "<<peak_neg/rms_neg<<" neg "<<saclst[isac]<<"\n";
      }
   }
   fout.close();

   return 0;
}
