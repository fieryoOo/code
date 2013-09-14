#define MAIN
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include "/home/tianye/code/Programs/head/mysac64.h"
#define fMax 2000
using namespace std;

void rmresponse_(int *n,double *dt,float *sei, double *freq, double *phase_res, double *amp_res, int *nf);

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


int main (int argc, char *argv[])
{
   if(argc != 4) {
      printf("Usage: do_cor [in_sac_csta] [in_sac_osta] [out_sac]\n");
      exit(-1);
   }
  float *sigc, *sigo;
  SAC_HD sdc, sdo;
  double pi = 4*atan(1.0);

  if( (read_sac(argv[1], &sigc, &sdc))==NULL ) {
      fprintf(stderr, "*** Warning: cannot read sac file %s ***", argv[1]);
      exit(0);
  }
  if( (read_sac(argv[2], &sigo, &sdo))==NULL ) {
      fprintf(stderr, "*** Warning: cannot read sac file %s ***", argv[2]);
      exit(0);
  }
  if(sdc.delta!=sdo.delta || sdc.npts!=sdc.npts) {
     cout<<"header file mismach!"<<endl;
     exit(0);
  }
  
  int i, bc, bo, dnum = 0;
  int lag = 500, nhalf = (int)floor(lag/sdc.delta);
  float cor[2*nhalf+1];
  for(i=0; i<nhalf; i++) {
     cor[i] = 0.;
     for(bc=nhalf-i,bo=0; bc<sdc.npts; bo++,bc++) cor[i] += sigc[bc] * sigo[bo];
     cor[i] /= sdc.npts-nhalf+i;
  }
  for(i=nhalf; i<=2*nhalf; i++) {
     cor[i] = 0.;
     for(bc=0,bo=i-nhalf; bo<sdo.npts; bo++,bc++) cor[i] += sigc[bc] * sigo[bo];
     cor[i] /= sdo.npts+nhalf-i;
  }

  strcpy(sdo.kevnm, sdc.kstnm);
  sdo.evla =  sdc.stla;
  sdo.evlo =  sdc.stlo;
  sdo.b = -nhalf*sdo.delta;
  sdo.npts = nhalf*2+1;
  sdo.user0 = dnum;
  sdo.nzjday = 1;
  write_sac (argv[3], cor, &sdo );

  return 1;  
}
