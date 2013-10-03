#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
using namespace std;

#define NARR 1000000

//--------Clock, swap & median----------//
uint64_t ClockGetTime()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
}

int main(int argc, char *argcv[])
{
   int i=10, nthreads = 10, niter = 15;
   float *a = new float[NARR], sum;
   double ts = ClockGetTime();
   //omp_set_dynamic(0);
   //omp_set_num_threads(nthreads);
   #pragma omp parallel for lastprivate(i) shared(a)
   //#pragma omp parallel
   for(i=0; i<niter; i++){
      int j, tid = omp_get_thread_num();
      if( tid == 0 ) {
	 nthreads = omp_get_num_threads();
	 printf("%d threads created.\n", nthreads);
      }
      fprintf(stdout, "Thread %d started at %f ms.  i = %d\n", tid, (ClockGetTime()-ts)/1.e3, i);
      for(j=0;j<NARR;j++){
	 //#pragma omp critical
	 #pragma omp atomic
	 a[j] += 1;
      }
      fprintf(stdout, "Thread %d ended at %f ms.  i = %d\n", tid, (ClockGetTime()-ts)/1.e3, i);
   }
   //#pragma omp barrier
   cout<<setprecision(15)<<"All threads done at "<<(ClockGetTime()-ts)/1.e3<<"ms.  i = "<<i<<endl;
   for(i=0;i<NARR;i++) if( a[i]!=niter ) cout<<i<<": "<<a[i]<<endl;
   delete [] a; a=NULL;
   return 0;
}
