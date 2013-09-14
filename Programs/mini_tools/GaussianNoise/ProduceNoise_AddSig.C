#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <time.h>
#include <sys/time.h>
#include "/home/tianye/code/Programs/head/mysac64.h"
using namespace std;

SAC_HD *read_shd (char *fname) {
   FILE *fsac;
   SAC_HD *SHD = new SAC_HD;
   if((fsac = fopen(fname, "r"))==NULL) return NULL;
   if ( !SHD ) SHD = &SAC_HEADER;
   fread(SHD,sizeof(SAC_HD),1,fsac);
   fclose (fsac);
   return SHD;
}

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

uint64_t ClockGetTime()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
}

void AddSig( float *sig, SAC_HD *shd, float A0, float per, float shift, float hw ) {
   double pi = 4*atan(1.);
   float alpha = 0.5/hw/hw, deltat;

   for(int i=0; i<shd->npts; i++) {
      deltat = i*shd->delta-shift;
      sig[i] += A0 * cos(2*pi/per*deltat) * exp(-alpha*pow(deltat,2));
   }
}

int main (int argc, char *argv[])
{
   if(argc != 7) {
      printf("Usage: Produce_noise [in_sac] [out_sac] [sig_per (-1 for noise)] [Tshift] [A0] [sig half hwidth]\n");
      exit(-1);
   }

   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::normal_distribution<> normal(0,1);
   boost::variate_generator< RNGType, boost::normal_distribution<> > Gauss(gener, normal);

   SAC_HD *shd = read_shd(argv[1]);
   if( shd == NULL ) {
      cout<<"Cannot access the header of "<<argv[1]<<endl;
      exit(0);
   }

   float per = atof(argv[3]), shift = atof(argv[4]), A0 = atof(argv[5]), hw = atof(argv[6]);
   if( (per<0&&per!=-1) || shift<0 || A0<0 || hw<0 ) {
      cout<<"Input parameters cannot be negative!"<<endl;
      exit(-1);
   }

   float *sig = new float[shd->npts];
   if( per == -1 ) {
      for(int i=0;i<shd->npts;i++) sig[i] = A0*Gauss();
      shd->user0 = -1;
   }
   else if (per <= 0.) {
      cout<<"non-positive per input!"<<endl;
      exit(0);
   }
   else {
      read_sac( argv[1], &sig, shd);
      AddSig( sig, shd, A0, per, shift, hw );
   }

   write_sac(argv[2], sig, shd);

   delete shd, sig;

   return 1;  
}
