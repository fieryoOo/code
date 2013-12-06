// multiple types of filters, applied to sac files
// bandpass: all fi > 0
// lowpass: f2=-1
// gaussian: f1=-1, f2=cper, f3=alpha, f4=-1

#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include "mysac64.h"
//#include <pthread.h>

#define PI 3.14159265358979323846

//extern pthread_mutex_t fftlock;

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD) {
   FILE *fsac;
   if((fsac = fopen(fname, "rb"))==NULL) return NULL;
   //pthread_mutex_lock(&fiolock);
   //if ( !SHD ) SHD = &SAC_HEADER;
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

void write_sac (char *fname, float *sig, SAC_HD *SHD) {
   FILE *fsac;
   if( (fsac = fopen(fname, "wb"))==NULL ) {
      std::cerr<<"ERROR(write_sac): Cannot open file "<<fname<<std::endl;
      return;
   }

   //if ( !SHD ) SHD = &SAC_HEADER;
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
    else if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }

//memset(SHD, 0, sizeof(SAC_HD));
	// pthread_mutex_lock(&fiolock);
         fwrite(SHD,sizeof(SAC_HD),1,fsac);
         fwrite(sig,sizeof(float),SHD->npts,fsac);
	// pthread_mutex_unlock(&fiolock);

        fclose (fsac);
}

void FFTW_F(fftw_plan plan, fftw_complex *out, int ns, float *seis, int n) {
   fftw_execute(plan);
   //pthread_mutex_lock(&fftlock);
   fftw_destroy_plan(plan);
   //pthread_mutex_unlock(&fftlock);
   int k;
   for(k=0; k<n; k++) seis[k] = out[k][0];
//for(k=0; k<n; k+=1000) if(seis[k] != 0.) printf("%d %f\n", k, seis[k]);
}

void FFTW_B(int type, float *seis, int n, fftw_complex **in, fftw_complex **out, int *nso, fftw_plan *planF, int Fflag) {
   int ns = (int)(log((double)n)/log(2.))+1;
   if(ns<13) ns = 13;
   ns = (int)pow(2,ns); *nso = ns;
   *in = (fftw_complex *) fftw_malloc ( ns * sizeof(fftw_complex) );//fftw_alloc_complex(ns);
   *out = (fftw_complex *) fftw_malloc ( ns * sizeof(fftw_complex) );//fftw_alloc_complex(ns);

   //measure plan using the allocated in/out blocks
   //pthread_mutex_lock(&fftlock);
   fftw_plan plan = fftw_plan_dft_1d (ns, *in, *out, FFTW_BACKWARD, type); //FFTW_ESTIMATE / FFTW_MEASURE
   if( Fflag == 1 ) *planF = fftw_plan_dft_1d (ns, *out, *in, FFTW_FORWARD, type);
   if( plan==NULL || (Fflag==1 && *planF==NULL) ) {
      fprintf(stderr,"Error(FFTW_B): fftw_plan creation failed!!\n");
      //pthread_mutex_unlock(&fftlock); exit(0);
   }
   //pthread_mutex_unlock(&fftlock);
   //initialize input array and excute
   memset(*in, 0, ns*sizeof(fftw_complex));
   int k;
   for(k=1; k<n; k++) (*in)[k][0] = seis[k];
   fftw_execute(plan);
   //cleanup
   //pthread_mutex_lock(&fftlock);
   fftw_destroy_plan(plan);
   //pthread_mutex_unlock(&fftlock);
   //if( Fflag==0 ) fftw_free(*in);

   //kill half spectrum and correct ends
   int nk = ns/2+1;
   for(k=nk;k<ns;k++) {
      (*out)[k][0] = 0.;
      (*out)[k][1] = 0.;
   }
   (*out)[0][0] /= 2.; (*out)[0][1] /= 2.;
   (*out)[nk-1][1] = 0.;
}

void TaperL( double f3, double f4, double dom, int nk, fftw_complex *sf ) {
   double f, ss;

   int i = (int)ceil(f3/dom);
   for(f=i*dom; f<f4; i++, f+=dom) {
      ss = ( 1. + cos(PI*(f3-f)/(f4-f3)) ) / 2.;
      sf[i][0] *= ss;
      sf[i][1] *= ss;
   }
   for(;i<nk;i++) {
      sf[i][0] = 0.;
      sf[i][1] = 0.;
   }

   return;
}

void TaperB( double f1, double f2, double f3, double f4, double dom, int nk, fftw_complex *sf ) {
   int i;
   double f, ss;

   for(i=0, f=0.; f<f1; i++, f+=dom) {
      sf[i][0] = 0.;
      sf[i][1] = 0.;
   }
   for(; f<f2; i++, f+=dom) {
      ss = ( 1. - cos(PI*(f1-f)/(f2-f1)) ) / 2.;
      sf[i][0] *= ss;
      sf[i][1] *= ss;
   }
   i = (int)ceil(f3/dom);
   for(f=i*dom; f<f4; i++, f+=dom) {
      ss = ( 1. + cos(PI*(f3-f)/(f4-f3)) ) / 2.;
      sf[i][0] *= ss;
      sf[i][1] *= ss;
   }
   for(;i<nk;i++) {
      sf[i][0] = 0.;
      sf[i][1] = 0.;
   }

   return;
}

void TaperGaussian( double fcenter, double fhlen, double dom, int nk, fftw_complex *sf ) {
   int i;
   float gauamp;
   double f, fstart, fend, fmax = (nk-1)*dom;
   // define effective window for given gaussian halflength
   f = fhlen * 4.;
   fstart = fcenter - f;
   fend = fcenter + f; 
   if( fend > fmax ) fend = fmax;
   //std::cerr<<"fcenter "<<fcenter<<"  dom "<<dom<<"  fstart "<<fstart<<"  fend "<<fend<<"  fmax "<<fmax<<std::endl;
   // cut high frequecies
   if( fstart > 0. ) {
      for(i=0, f=0.; f<fstart; i++, f+=dom) {
	 sf[i][0] = 0.;
	 sf[i][1] = 0.;
      }
      f = i * dom; // correct for round-off error
   }
   else { f = 0.; i = 0; }
   // apply taper
   float alpha = -0.5/(fhlen*fhlen);
   for(; f<fend-1.e-10; i++, f+=dom) {
      gauamp = f - fcenter;
      gauamp = exp( alpha * gauamp * gauamp );
      sf[i][0] *= gauamp;
      sf[i][1] *= gauamp;
   }
   // cut low frequencies
   if( fend < fmax ) {
      f = i * dom; // again, correct for round-off
      for(; i<nk; i++) {
         sf[i][0] = 0.;
         sf[i][1] = 0.;
      }
   }
   
}

void Filter (double f1, double f2, double f3, double f4, double dt, int n, float *seis_in, float *seis_out) {
   if(f4 > 0.5/dt) {
      fprintf(stdout, "\n*** Warning(Filter): filter band out of range! ***\n");
      f4 = 0.49999/dt;
   }
   fftw_plan planF = NULL;
   //backward FFT: s ==> sf
   int ns;
   fftw_complex *s, *sf;
   FFTW_B(FFTW_ESTIMATE, seis_in, n, &s, &sf, &ns, &planF, 1);
   //make tapering
   int nk = ns/2+1;
   double dom = 1./dt/ns;
   if( (f1==-1. || f2==-1.) && (f3>0. && f4>0.) ) TaperL( f3, f4, dom, nk, sf );
   else if( f1>=0 && f2>0 && f3>0 && f4>0 ) TaperB( f1, f2, f3, f4, dom, nk, sf );
   else if( f1==-1. && f4==-1. ) TaperGaussian( f2, f3, dom, nk, sf );
   else {
      std::cerr<<"Unknow filter type for parameters = "<<f1<<" "<<f2<<" "<<f3<<" "<<f4<<std::endl;
      exit(-1);
   }

   //forward FFT: sf ==> s
   FFTW_F(planF, s, ns, seis_out, n);
   fftw_free(s); fftw_free(sf);

   //forming final result
   int k;
   float ftmp = 2./ns;
   for(k=0; k<n; k++) {
      if( seis_in[k]==0 ) seis_out[k] = 0.;
      else seis_out[k] *= ftmp;
   }

   return;
}


int main(int argc, char *argv[])
{
   if(argc != 6) {
      std::cerr<<"Usage: "<<argv[0]<<" [sac file] [f1 (-1=gaussian)] [f2 (-1=lowpass fcenter=gaussian)] [f3 (frec halflength=gaussian)] [f4 (-1=gaussian)]"<<std::endl;
      exit(-1);
   }

   SAC_HD shd;
   float *sig=NULL;
   if (read_sac(argv[1], &sig, &shd)==NULL) return 0;

   double f1 = atof(argv[2]), f2 = atof(argv[3]), f3 = atof(argv[4]), f4 = atof(argv[5]);
   Filter(f1, f2, f3, f4, (double)shd.delta, shd.npts, sig, sig);

   char outname[300];
   sprintf(outname, "%s_ft", argv[1]);
   write_sac(outname, sig, &shd);

   return 0;
}
