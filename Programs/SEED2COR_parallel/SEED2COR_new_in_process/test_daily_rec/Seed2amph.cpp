#include "SysTools.h"
#include "DailyRec.h"
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <vector>
#include <string>
#include "MyOMP.h"

/*
void SetName(int ne, int ns) {
   sprintf(sdb->rec[ne][ns].fname,"%s/%s/%s.%s.%s.SAC", sdb->mo[imonth].name, sdb->ev[ne].name, sdb->ev[ne].name, sdb->st[ns].name, ch);
   sprintf(sdb->rec[ne][ns].ft_fname,"%s/%s/ft_%s.%s.%s.SAC", sdb->mo[imonth].name, sdb->ev[ne].name, sdb->ev[ne].name, sdb->st[ns].name, ch);
   sprintf(sdb->rec[ne][ns].chan,"%s", ch );
}
*/



int  main ( int argc, char *argv[] ) 
{
   if( argc != 2) {
      std::cerr<<"Usage: "<<argv[0]<<" [input parameter file]"<<std::endl;
      return -1;
   }

   /* store old SAC and RESP files if there's any */
   //bool oldfiles = false;
   MKDir("old_sac_files");
   //wMove(".", "*.SAC", "old_sac_files", 0, &nfmvd); if( nfmvd > 0 ) oldfiles = true;
   //wMove(".", "RESP.*", "old_sac_files", 0, &nfmvd); if( nfmvd > 0 ) oldfiles = true;
   std::vector<std::string> filelst;
   bool oldsac = wMove(".", "*.SAC", "old_sac_files", filelst);
   bool oldresp = wMove(".", "RESP.*", "old_sac_files", filelst);


   //int ithread; 
   /* create the DailyRec object, load in parameters from the input file */
   DailyRec dailyrec(argv[1]);

   /* set/check parameters and extract osac */
   //dailyrec.Set("sps 20");
   if( ! dailyrec.CheckPreExtract() ) return -1;
   std::vector<std::string> stalst = {"109C", "A04A", "ACSO", "ADO"};
   #pragma omp parallel for
   for(int i=0; i<stalst.size(); i++) {
      // copy params from dailyrec
      DailyRec DRtmp(dailyrec);
      // set staname
      std::string stmp = "staname " + stalst.at(i);
      DRtmp.Set(stmp.c_str());
      // set fosac
      stmp = "fosac " + stalst.at(i) + ".sac";
      DRtmp.Set(stmp.c_str());
      // set ffsac
      stmp = "ffsac " + stalst.at(i) + "_ft.sac";
      DRtmp.Set(stmp.c_str());

      // extract sac, donot skip | write to disc
      if( ! DRtmp.ExtractSac(0, true) ) continue;
      //DRtmp.RmRESP(true);
   }
  
   if( ! dailyrec.CheckPreRmRESP() ) return -1;

   if( ! dailyrec.CheckPreTSNorm() ) return -1;

   /* fetch back old files and remove temporary dir */
   /*
   if( oldfiles ) {
      wMove("old_sac_files", "*.SAC", ".", 0, &nfmvd);
      wMove("old_sac_files", "RESP.*", ".", 0, &nfmvd);
   }
   */
   if( oldsac ) wMove("old_sac_files", "*.SAC", ".", filelst);
   if( oldresp ) wMove("old_sac_files", "RESP.*", ".", filelst);
   fRemove("old_sac_files");

   return 0;
}
