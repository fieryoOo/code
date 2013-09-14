#define MAIN
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include "/home/tianye/code/Programs/head/mysac64.h"
using namespace std;

#define PI 3.1415926536

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD) {
   FILE *fsac;
   if((fsac = fopen(fname, "r"))==NULL) return NULL;
   if ( !SHD ) SHD = &SAC_HEADER;
   fread(SHD,sizeof(SAC_HD),1,fsac);
   *sig = (float *) malloc (SHD->npts * sizeof(float));
   fread(*sig,sizeof(float),SHD->npts,fsac);
   fclose (fsac);

   /*-------------  calcule de t0  ----------------*/
   {
        int eh, em ,i;
        float fes;
        char koo[9];

        for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
        koo[8] = 0;

        SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
         SHD->nzsec + SHD->nzmsec*.001;

        sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

        SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}
   return SHD;
}

void write_sac (char *fname, float *sig, SAC_HD *SHD) {
   FILE *fsac;
   fsac = fopen(fname, "wb");
   if ( !SHD ) SHD = &SAC_HEADER;
   SHD->iftype = (int)ITIME;
   SHD->leven = (int)TRUE;
   SHD->lovrok = (int)TRUE;
   SHD->internal4 = 6L;

  /*+++++++++++++++++++++++++++++++++++++++++*/
   SHD->depmin = sig[0];
   SHD->depmax = sig[0];
   int i;
   for ( i = 0; i < SHD->npts ; i++ )
   {
    if ( SHD->depmin > sig[i] ) SHD->depmin = sig[i];
    if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }

         fwrite(SHD,sizeof(SAC_HD),1,fsac);

         fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);


        fclose (fsac);
}

int nint (float a) {
   return (int)floor(a+0.5);
}

int main (int argc, char *argv[])
{
   if(argc != 8) {
      printf("Usage: %s [in_sacf] [out_name] [tll] [tlh] [thl] [thh] [sym?]\n", argv[0]);
      exit(-1);
   }

   float *sig;
   SAC_HD shd;
   if( (read_sac(argv[1], &sig, &shd))==NULL ) {
      fprintf(stderr, "*** Warning: cannot read sac file %s ***\n", argv[1]);
      exit(0);
  }
  float tll = atof(argv[3]), tlh = atof(argv[4]), thl = atof(argv[5]), thh = atof(argv[6]);
  int i, nb, ne;
  float theta, step;
  int fsym = atoi(argv[7]);
  if( fsym ) {
     nb = nint((-thh-shd.b)/shd.delta); ne = nint((-thl-shd.b)/shd.delta);
     for(i=0;i<nb;i++) sig[i] = 0.;
     for(i=nb,theta=0.,step=PI/(ne-nb);i<ne;i++,theta+=step) sig[i] *= 0.5*(1.-cos(theta));
     nb = nint((-tlh-shd.b)/shd.delta); ne = nint((-tll-shd.b)/shd.delta);
     for(i=nb,theta=-PI,step=PI/(ne-nb);i<ne;i++,theta+=step) sig[i] *= 0.5*(1.-cos(theta));
  }
  else ne = 0;
  nb = nint((tll-shd.b)/shd.delta);
  for(i=ne;i<nb;i++) sig[i] = 0.;
  ne = nint((tlh-shd.b)/shd.delta);
  for(i=nb,theta=0.,step=PI/(ne-nb);i<ne;i++,theta+=step) sig[i] *= 0.5*(1.-cos(theta));
  nb = nint((thl-shd.b)/shd.delta); ne = nint((thh-shd.b)/shd.delta);
  for(i=nb,theta=-PI,step=PI/(ne-nb);i<ne;i++,theta+=step) sig[i] *= 0.5*(1.-cos(theta));
  for(i=ne;i<shd.npts;i++) sig[i] = 0.;

  write_sac(argv[2], sig, &shd);
  return 1;
}
