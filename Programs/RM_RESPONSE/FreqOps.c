#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>
//#include <pthread.h>

#define PI 3.14159265358979323846

//extern pthread_mutex_t fftlock;

void FFTW_F(fftw_plan *plan, int type, fftw_complex *in, int ns, float *seis, int n) {
   fftw_complex *out = (fftw_complex *) calloc ( ns, sizeof(fftw_complex) );
   //pthread_mutex_lock(&fftlock);
   *plan = fftw_plan_dft_1d (ns, in, out, FFTW_FORWARD, type); //FFTW_ESTIMATE / FFTW_MEASURE
   //pthread_mutex_unlock(&fftlock);
   fftw_execute(*plan);
   //pthread_mutex_lock(&fftlock);
   fftw_destroy_plan(*plan);
   //pthread_mutex_unlock(&fftlock);
   for(int k=0; k<n; k++) seis[k] = out[k][0];
   free(out);
}

void FFTW_B(fftw_plan *plan, int type, float *seis, int n, fftw_complex **out, int *nso) {
   int ns = (int)(log((double)n)/log(2.))+1;
   if(ns<13) ns = 13;
   ns = (int)pow(2,ns); *nso = ns;
   fftw_complex *in = (fftw_complex *) calloc ( ns, sizeof(fftw_complex) );
   *out = (fftw_complex *) malloc ( ns * sizeof(fftw_complex) ); 

   int k;
   for(k=1; k<n; k++) in[k][0] = seis[k];

   //pthread_mutex_lock(&fftlock);
   *plan = fftw_plan_dft_1d (ns, in, *out, FFTW_BACKWARD, type); //FFTW_ESTIMATE / FFTW_MEASURE
   //pthread_mutex_unlock(&fftlock);
   fftw_execute(*plan);
   //pthread_mutex_lock(&fftlock);
   fftw_destroy_plan(*plan);
   //pthread_mutex_unlock(&fftlock);

   //kill half spectrum and correct ends
   int nk = ns/2+1;
   for(k=nk;k<ns;k++) {
      (*out)[k][0] = 0.;
      (*out)[k][1] = 0.;
   }
   (*out)[0][0] /= 2.; (*out)[0][1] /= 2.;
   (*out)[nk-1][1] = 0.;
   free(in);
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

void fDiv(double dom, int nsig, fftw_complex *sf, double *freq, double *amp, double *pha, int ntra) {
   int isig, itra = 1;
   double f, ampcur, phacur, realsig, imagsig;
   double sintmp, costmp;
   for(isig=(int)ceil(freq[0]/dom); isig<nsig; isig++) {
      f = isig*dom;
      while(f > freq[itra]) {
	 itra++;
	 if(itra >= ntra) break;
      }
      //interpolate to get current amp and pha
      sintmp = (f-freq[itra-1]) / (freq[itra]-freq[itra-1]);
      ampcur = amp[itra-1] + (amp[itra]-amp[itra-1]) * sintmp;
      phacur = pha[itra-1] + (pha[itra]-pha[itra-1]) * sintmp;
      //divide sf by (ampcur, phacur)
      realsig = sf[isig][0]; imagsig = sf[isig][1];
      sintmp = sin(phacur); costmp = cos(phacur);
      sf[isig][0] = (realsig*costmp - imagsig*sintmp) / ampcur;
      sf[isig][1] = (imagsig*costmp + realsig*sintmp) / ampcur;
   }

}

void FDivide (double f1, double f2, double f3, double f4, double dt, int n, float *seis_in, float *seis_out, double *freq, double *amp, double *pha, int nf) {
   if(f4 > 0.5/dt) {
      fprintf(stderr, "### Warning: filter band out of range! ###");
      return;
   }
   fftw_plan plan1;
   //backward FFT: s ==> sf
   int ns;
   fftw_complex *sf;
   FFTW_B(&plan1, FFTW_ESTIMATE, seis_in, n, &sf, &ns);
   //make tapering
   int nk = ns/2+1;
   double dom = 1./dt/ns;
   fDiv(dom, nk, sf, freq, amp, pha, nf);
   TaperB( f1, f2, f3, f4, dom, nk, sf );

   //forward FFT: sf ==> s
   fftw_plan plan2;
   FFTW_F(&plan2, FFTW_ESTIMATE, sf, ns, seis_out, n);
   free(sf);

   //forming final result
   int k;
   for(k=0; k<n; k++) {
      if( seis_in[k]==0 ) seis_out[k] = 0.;
      else seis_out[k] *= 2./ns;
   }

   return;
}

