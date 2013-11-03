#include <stdio.h>
#include <iostream>
#include <cmath>
#include "mysac64.h"

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD) {
   FILE *fsac;
   if((fsac = fopen(fname, "rb"))==NULL) return NULL;
   //pthread_mutex_lock(&fiolock);
   if ( !SHD ) SHD = &sac_null;
   fread(SHD,sizeof(SAC_HD),1,fsac);
   if( *sig == NULL) *sig = (float *) malloc (SHD->npts * sizeof(float));
   //else *sig = (float *) realloc (*sig, SHD->npts * sizeof(float));
   fread(*sig,sizeof(float),SHD->npts,fsac);
   fclose (fsac);
   //pthread_mutex_unlock(&fiolock);

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

int main(int argc, char *argv[])
{
   if( argc != 4 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [in sac file] [grv min] [grv max]"<<std::endl;
      exit(-1);
   }
   float grvmin=atof(argv[2]), grvmax=atof(argv[3]);
   if( grvmin<0. || grvmax>10. || grvmin>=grvmax ) {
      std::cerr<<"Incorrect input-group-vel range: "<<grvmin<<" - "<<grvmax<<std::endl;
      return -1;
   }

   // read in signal
   float *sig;
   SAC_HD shd;
   if (read_sac(argv[1], &sig, &shd)==NULL) {
      std::cerr<<"Connot read sac file "<<argv[1]<<std::endl;
      return -1;
   }
   if( shd.b!=0 ) {
      std::cerr<<"Sac file doesn't begin at 0. Modifications needed!"<<std::endl;
      return -1;
   }

   // define signal time window
   float tb = shd.dist/grvmax, te = shd.dist/grvmin;
   int i, idxb = (int)floor(tb/shd.delta+0.5), idxe = (int)floor(te/shd.delta+0.5) + 1;
   if( idxe > shd.npts ) idxe = shd.npts;
   // compute averaged amplitude in given window
   float amprms=0.;
   for(i=idxb; i<idxe; i++) amprms += sig[i] * sig[i];
   amprms = sqrt(amprms/(idxe-idxb-1));

   // define noise
   tb = shd.e-1000.; 
   if( tb < te+500 ) tb = te + 500.;
   te = tb + 500.;;
   if( te > shd.e - 100. ) te = shd.e - 100;
   if( (te-tb) < 100. ) {
      std::cerr<<"No enough trailing points for computing noise"<<std::endl;
      return -1;
   }
   // compute rms noise
   idxb = (int)floor(tb/shd.delta+0.5); idxe = (int)floor(te/shd.delta+0.5) + 1;
   float noiserms = 0.;
   for(i=idxb;i<idxe;i++) noiserms += sig[i] * sig[i];
   noiserms=sqrt(noiserms/(idxe-idxb-1));

   std::cout<<amprms<<" "<<noiserms<<" "<<amprms/noiserms<<std::endl;
   free(sig);
   return 0;
}
