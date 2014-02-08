#include "SysTools.h"
#include "DailyRec.h"
#include <cstdlib>
#include <cstring>

/*
void SetName(int ne, int ns) {
   sprintf(sdb->rec[ne][ns].fname,"%s/%s/%s.%s.%s.SAC", sdb->mo[imonth].name, sdb->ev[ne].name, sdb->ev[ne].name, sdb->st[ns].name, ch);
   sprintf(sdb->rec[ne][ns].ft_fname,"%s/%s/ft_%s.%s.%s.SAC", sdb->mo[imonth].name, sdb->ev[ne].name, sdb->ev[ne].name, sdb->st[ns].name, ch);
   sprintf(sdb->rec[ne][ns].chan,"%s", ch );
}
*/



int  main () 
{
   /* store old SAC and RESP files if there's any */
   bool oldfiles = false;
   MKDir("old_sac_files");
   int nfmvd;
   wMove(".", "*.SAC", "old_sac_files", 0, &nfmvd); if( nfmvd > 0 ) oldfiles = true;
   wMove(".", "RESP.*", "old_sac_files", 0, &nfmvd); if( nfmvd > 0 ) oldfiles = true;

   //int ithread; 
   /* create the DailyRec object */
   std::string rdsexe("/home/yeti4009/usr/bin/rdseed"), evrexe("/home/yeti4009/usr/bin/evalresp"), fseed("./OBS_US_NEW_2012.JAN.1.457078.seed");
   std::string staname("I05D"), chname("BHZ");
   std::string fosac("./temp.sac"), ffsac("./temp_ft.sac"), outdir(".");
   DailyRec dailyrec(rdsexe, evrexe, fseed, staname, chname, fosac, ffsac, outdir);
   dailyrec.ExtractSac();
  
   /* fetch back old files and remove temporary dir */
   if( oldfiles ) {
      wMove("old_sac_files", "*.SAC", ".", 0, &nfmvd);
      wMove("old_sac_files", "RESP.*", ".", 0, &nfmvd);
   }
   fRemove("old_sac_files");

   return 0;
}
