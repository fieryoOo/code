#include <stdio.h>
#include <stdlib.h>
#include "/home/tianye/code/Programs/head/mysac.h"

SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*--------------------------------------------------------------------------*/
/* function to read sac files given the name, fname. The function outputs the time signal to the pointer sig
, fills the header SHD, if the signal has fewer than nmax points */
{
  FILE *fsac;

  if((fsac = fopen(fname, "rb")) == NULL) {
    fprintf(stderr,"read_sac: Could not open %s\n", fname);
    return NULL;
  }

  if ( !fsac ) {
    /*fprintf(stderr,"file %s not find\n", fname);*/
    return NULL;
  }

//  if ( !SHD ) SHD = &SAC_HEADER;

  fread(SHD,sizeof(SAC_HD),1,fsac);

  if ( SHD->npts > nmax ) {
    fprintf(stderr,"ATTENTION !!! %s npts is limited to %d.\n", fname, (int)nmax);
    SHD->npts = nmax;
  }

  fread(sig,sizeof(float),(int)(SHD->npts),fsac);
  fclose (fsac);

/*-------------  calcule de t0  ----------------*/
   {
     int eh, em ,i;
     float fes;
     char koo[9];

     for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
     koo[8] = '\0';

     SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
     SHD->nzsec + SHD->nzmsec*.001;

     sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

     SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}

   return SHD;
}


void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*--------------------------------------------------------------------------*/
{
  FILE *fsac;
  int i;
  if((fsac = fopen(fname, "wb"))==NULL) {
    fprintf(stderr,"write_sac: Could not open %s to write\n", fname);
  }
  else {

    if ( !SHD ) {
//      SHD = &SAC_HEADER;
    }

    SHD->iftype = (int)ITIME;
    SHD->leven = (int)TRUE;
    SHD->lovrok = (int)TRUE;
    SHD->internal4 = 6L;
    SHD->depmin = sig[0];
    SHD->depmax = sig[0];

    for ( i = 0; i < SHD->npts ; i++ ) {
      if ( SHD->depmin > sig[i] ) {
        SHD->depmin = sig[i];
      }
      if ( SHD->depmax < sig[i] ) {
        SHD->depmax = sig[i];
      }
    }

    fwrite(SHD,sizeof(SAC_HD),1,fsac);
    fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);

    fclose (fsac);
  }
}

