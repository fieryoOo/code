#ifndef SACREC_H
#define SACREC_H

#include "mysac64.h"
#include <cstdio>
#include <string>
#include <iostream>
#include <tr1/memory>

class SacRec {
private:
   std::string fname;
   SAC_HD shd;
   std::shared_ptr<float> sig;
public:
   /* read sac header */
   bool read_shd ();
   /* read sac file, signal memory is allocated on heap */
   SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD);
   /* write sac file */
   void write_sac (const char *fname, float *sig, SAC_HD *SHD);

   void UpdateTime();

};

#endif
