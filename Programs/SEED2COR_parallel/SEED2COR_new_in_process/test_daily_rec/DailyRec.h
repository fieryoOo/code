#ifndef DAILYREC_H
#define DAILYREC_H

#include "mysac64.h"
#include <iostream>
#include <cstring>

#ifndef PathLen
#define PathLen 200
#endif

class DailyRec {
private:
   std::string rdsexe, evrexe;
   std::string fseed, fresp, fosac, ffsac;
   int fskipesac;
   int npts;

public:
   /* read in parameters through constructor */
   DailyRec( std::string rdsexein, std::string evrexein, std::string fseedin, std::string fosacin, std::string ffsacin )
      : rdsexe(rdsexein), evrexe(evrexein), fseed(fseedin), fosac(fosacin), ffsac(ffsacin) {}
   // list of record objects with information of event(fseed), station(lon & lat) and singlestation-files(f_original_sac, f_ft_sac, f_am, f_ph)
   int CheckExistence(int ithread);
   void ExtractSac();
   void RmRESP();
   void TempSpecNorm ();
};


#endif
