#include <stdlib.h>
#include <fftw3.h>
//#include "/home/tianye/src/local/include/fftw3.h"
#include <math.h>
#include <string.h>
#include <pthread.h>

extern pthread_mutex_t fftlock;

void FFTW_F(fftw_plan *plan, int type, fftw_complex *in, int ns, float *seis, int n);

void IFFT(int nlen, float *amp, float *pha, float *seis_out, int *nsig) {
   int i, ns = (nlen-1)*2;
   fftw_plan plan3;
   fftw_complex *in;// *out;

   in = (fftw_complex *) calloc ( ns, sizeof(fftw_complex) );
   for(i=0;i<nlen;i++) {
      in[i][0] = amp[i] * cos(pha[i]);
      in[i][1] = -amp[i] * sin(pha[i]);
   }
   FFTW_F(&plan3, FFTW_MEASURE, in, ns, seis_out, ns);

   free(in);
   *nsig = ns;

   return;
}

