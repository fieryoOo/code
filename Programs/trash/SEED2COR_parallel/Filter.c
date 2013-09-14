#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

#define PI 3.14159265358979323846

extern pthread_mutex_t fftlock;

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

void Filter (double f1, double f2, double f3, double f4, double dt, int n, float *seis_in, float *seis_out) {
   if(f4 > 0.5/dt) {
      fprintf(stderr, "### Warning: filter band out of range! ###");
      return;
   }
   //determin the power of FFT
   int ns = (int)pow(2., (int)(log((double)n)/log(2.))+1);
   double dom = 1./dt/ns;
   int k;

   //arrays and plans for FFT
   fftw_plan plan1, plan2;
   fftw_complex *s, *sf;
   s = (fftw_complex *) calloc ( ns, sizeof(fftw_complex) );
   sf = (fftw_complex *) malloc ( ns * sizeof(fftw_complex) );

   for(k=0;k<n;k++) s[k][0] = seis_in[k];

   //backward FFT: s ==> sf
   pthread_mutex_lock(&fftlock);
   plan1 = fftw_plan_dft_1d (ns, s, sf, FFTW_BACKWARD, FFTW_ESTIMATE); //FFTW_ESTIMATE / FFTW_MEASURE
   pthread_mutex_unlock(&fftlock);
   fftw_execute(plan1);
   pthread_mutex_lock(&fftlock);
   fftw_destroy_plan(plan1);
   pthread_mutex_unlock(&fftlock);

   //kill half spectrum and correct ends
   int nk = ns/2+1;
   for(k=nk;k<ns;k++) {
      sf[k][0] = 0.;
      sf[k][1] = 0.;
   }
   sf[0][0] /= 2.; sf[0][1] /= 2.;
   sf[nk-1][1] = 0.;

   //make tapering
   if( f2 == -1. ) TaperL( f3, f4, dom, nk, sf );
   else TaperB( f1, f2, f3, f4, dom, nk, sf );

   //forward FFT: sf ==> s
   pthread_mutex_lock(&fftlock);
   plan2 = fftw_plan_dft_1d (ns, sf, s, FFTW_FORWARD, FFTW_ESTIMATE); //FFTW_ESTIMATE / FFTW_MEASURE
   pthread_mutex_unlock(&fftlock);
   fftw_execute(plan2);
   pthread_mutex_lock(&fftlock);
   fftw_destroy_plan(plan2);
   pthread_mutex_unlock(&fftlock);

   //forming final result
   for(k=0; k<n; k++)
      if( seis_in[k]==0 ) seis_out[k] = 0.;
      else seis_out[k] = s[k][0]*2./ns;

   free(s); free(sf);

   return;
}
