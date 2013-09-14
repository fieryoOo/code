#include <stdlib.h>
#include <fftw3.h>
//#include "/home/tianye/src/local/include/fftw3.h"
#include <math.h>
#include <string.h>
#include <pthread.h>

extern pthread_mutex_t fftlock;

void IFFT(int nlen, float *amp, float *pha, float *seis_out, int *nsig) {
   int i, ns = (nlen-1)*2;
   fftw_plan plan3;
   //fftw_complex in[ns], out[ns];
   fftw_complex *in, *out;

   in = (fftw_complex *) calloc ( ns, sizeof(fftw_complex) );
   out = (fftw_complex *) malloc ( ns * sizeof(fftw_complex) );
   for(i=0;i<nlen;i++) {
      in[i][0] = amp[i] * cos(pha[i]);
      in[i][1] = -amp[i] * sin(pha[i]);
   }

   pthread_mutex_lock(&fftlock);
   plan3 = fftw_plan_dft_1d (ns, in, out, FFTW_FORWARD, FFTW_MEASURE); //FFTW_ESTIMATE / FFTW_MEASURE
   pthread_mutex_unlock(&fftlock);
   fftw_execute(plan3);
   pthread_mutex_lock(&fftlock);
   fftw_destroy_plan(plan3);
   pthread_mutex_unlock(&fftlock);
   for(i=0;i<ns;i++) seis_out[i] = out[i][0];

   free(in); free(out);
   *nsig = ns;

   return;
}

int Whiten (double f1, double f2, double f3, double f4, int npow, double dt, int n, float hlen, float seis_in[], float seissm[], float seis_out[], float seis_outamp[], float seis_outph[], int *ns, double *dom, int *flag_whiten) {

}
