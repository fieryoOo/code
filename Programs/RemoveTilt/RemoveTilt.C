#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "mysac64.h"
using namespace std;

extern "C" {
void fft_(double *dt, int *n, float seis_in[], float seis_outamp[], float seis_outph[],int *nk,double *dom);
void rmresponse_(double *f1, double *f2, double *f3, double *f4,
   float seis[], double *freq, double *phase_res, double *amp_res);
//     int *n, double *dt, float seis[], double *freq, double *phase_res, double *amp_res);
}

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD) {
   FILE *fsac;
   if((fsac = fopen(fname, "r"))==NULL) return NULL;
   if ( !SHD ) SHD = &SAC_HEADER;
   fread(SHD,sizeof(SAC_HD),1,fsac);
   //*sig = (float *) malloc (SHD->npts * sizeof(float));
   *sig = new float[SHD->npts];
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
   if( argc != 8) {
      cout<<"Usage: "<<argv[0]<<" [SAC_BHZ] [SAC_BDH] [f1] [f2] [f3] [f4] [out_name]"<<endl;
      exit(-1);
   }
   
   double f1 = atof(argv[3]);
   double f2 = atof(argv[4]);
   double f3 = atof(argv[5]);
   double f4 = atof(argv[6]);
   if( f1<0 || f1>=f2 || f2>=f3 || f3>=f4 ) {
      cout<<"Incorrect corner frequencies"<<endl;
      exit(0);
   }

   float *sigZ, *sigD;
   SAC_HD shdZ, shdD;
   if( read_sac(argv[1], &sigZ, &shdZ) == NULL ) {
      cout<<"Cannot access file "<<argv[1]<<endl;
      exit(0);
   }
   if( read_sac(argv[2], &sigD, &shdD) == NULL ) {
      cout<<"Cannot access file "<<argv[2]<<endl;
      exit(0);
   }

   int nkZ, nkD, npow = 1;
   int nsp = (int)pow(2., floor(log(shdZ.npts)/log(2.)))+10;
   float sigZam[nsp], sigZph[nsp];
   float sigDam[nsp], sigDph[nsp];
   double dt = shdZ.delta, domZ, domD;
   if(dt!=shdD.delta) {
      cout<<"Incompatible sampling rate"<<endl;
      exit(0);
   }
   fft_(&dt,&shdZ.npts,sigZ,sigZam,sigZph,&nkZ,&domZ);
   fft_(&dt,&shdD.npts,sigD,sigDam,sigDph,&nkD,&domD);
   if( domZ!=domD || nkZ!=nkD ) {
      cout<<"Incompatible am/ph files"<<endl;
      exit(0);
   }

   int i, j, ii;
   int hlen = 50, ib = (int)floor(f1/domZ);
   if(ib<hlen+1) { ib=hlen; cout<<"Warning: Starting frequency too small!"<<endl; }
   double xx=0.,yy=0.,xy_r=0.,xy_i=0.,xy1,xy2;
   double temp_ph, temp_am, ftmp, fcurrent = ib*domZ, pi = 4.0*atan(1.0);
   double freq[nkZ], phres[nkZ], amres[nkZ];
   FILE *fout;
   char nametmp[100];
   sprintf(nametmp, "%s_TF", argv[7]);
   fout = fopen(nametmp, "w");
   for(i=ib-hlen-1; i<ib+hlen; i++) {
      xx += pow(sigDam[i],2.);
      yy += pow(sigZam[i],2.);
      xy1 = sigDam[i]*sigZam[i];
      xy_r += xy1*cos(sigDph[i]-sigZph[i]);
      xy_i += xy1*sin(sigDph[i]-sigZph[i]);
   }
   for(i=ib,ftmp=f4+domZ; i<nkZ-hlen&&fcurrent<=ftmp; i++, fcurrent+=domZ) {
      xx += pow(sigDam[i+hlen],2.) - pow(sigDam[i-hlen-1],2.);
      yy += pow(sigZam[i+hlen],2.) - pow(sigZam[i-hlen-1],2.);
      xy1 = sigDam[i+hlen]*sigZam[i+hlen]; xy2 = sigDam[i-hlen-1]*sigZam[i-hlen-1];
      xy_r += xy1 * cos(sigDph[i+hlen]-sigZph[i+hlen]) - xy2 * cos(sigDph[i-hlen-1]-sigZph[i-hlen-1]);
      xy_i += xy1 * sin(sigDph[i+hlen]-sigZph[i+hlen]) - xy2 * sin(sigDph[i-hlen-1]-sigZph[i-hlen-1]);
      temp_ph=atan(xy_i/xy_r);
      if(xy_r<0) temp_ph+=pi;
      if(temp_ph>pi) temp_ph-=2*pi;
      temp_am=sqrt(xy_r*xy_r+xy_i*xy_i);
      ii = i-ib;
      freq[ii] = fcurrent;
      phres[ii] = temp_ph;
      amres[ii] = temp_am/xx;
      fprintf(fout, "%g %g %g %g\n", fcurrent, temp_am/sqrt(xx*yy), temp_am/xx, temp_ph);
   }
   fclose(fout);

   dt = shdD.delta;
   rmresponse_(&f1,&f2,&f3,&f4,sigD,freq,phres,amres);
   //rmresponse_(&f1,&f2,&f3,&f4,&(shdD.npts),&dt,sigD,freq,phres,amres);
   //write_sac ("PRED_tmp.SAC", sigD, &shdD);
   for(i=0;i<shdZ.npts;i++) sigZ[i] -= sigD[i];
   write_sac (argv[7], sigZ, &shdZ );

   delete [] sigZ;
   delete [] sigD;
   return 1;
}
