#ifndef DAILYREC_H
#define DAILYREC_H

#include "SacRec.h"
#include <iostream>
#include <sstream>
#include <cstring>

#ifndef PathLen
#define PathLen 200
#endif

class DailyRec {
private:
   class DRimpl;
   std::unique_ptr<DRimpl> pimpl;

/*
   std::string rdsexe, evrexe;
   std::string fseed, fresp, fosac, ffsac;
   std::string staname, chname, outdir;
   std::ostringstream reports;
   char *resp_fname;
   float lon, lat;
   int year, month, day;
   int fskipesac;
   int npts;
   float t0, delta;
*/

public:
   /* read in parameters through constructor */
   DailyRec( std::string rdsexein, std::string evrexein, std::string fseedin, std::string stanamein, std::string chnamein, std::string fosacin, std::string ffsacin, std::string outdirin );
//      : rdsexe(rdsexein), evrexe(evrexein), fseed(fseedin), staname(stanamein), chname(chnamein), fosac(fosacin), ffsac(ffsacin), outdir(outdirin), pimpl(new DRimpl()) {}
   ~DailyRec();
   // list of record objects with information of event(fseed), station(lon & lat) and singlestation-files(f_original_sac, f_ft_sac, f_am, f_ph)
   //bool CheckExistence();
   void ExtractSac();
   //void RmRESP();
   //void TempSpecNorm ();
};


#endif
