#include <stdio.h>
#include <string.h>
#define MAIN
#include "mysac.h"

SAC_HD shd;
/* function prototyps     */
void write_sac (char *fname, float *sig, SAC_HD *SHD);

int wsac_(char *fname, float *elat, float *elon, char *cod, float *slat,
          float *slon, float *dt, char *chn, char *net, int *n, float *seis)
{
int i;
  shd = sac_null;
  shd.delta = *dt;
  shd.npts = *n;
  shd.b = 0;
  shd.e = shd.delta*(shd.npts-1);
  shd.evdp = 0.0;
  shd.evel = 0.0;
  shd.evla = *elat;
  shd.evlo = *elon;
  strncpy(shd.kstnm,cod,8);
  shd.stdp = 0.0;
  shd.stel = 0.0;
  shd.stla = *slat;
  shd.stlo = *slon;
  shd.nzyear = 2000;
  shd.nzjday = 1;
  shd.nzhour = 0;
  shd.nzmin = 0;
  shd.nzsec = 0;
  shd.nzmsec = 0;
  strncpy(shd.kcmpnm,chn,8);
  strncpy(shd.knetwk,net,8);
  for(i=0; i < 256; i++) if(fname[i] == ' '){fname[i] = 0; break;};
  write_sac (fname, seis, &shd);
  return 0;
}
