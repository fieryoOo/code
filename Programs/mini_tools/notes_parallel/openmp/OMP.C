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
   int i=0, nthreads = 20;
   float *a = new float[nthreads], *b = new float[NARR];
   double ts = ClockGetTime();
   //omp_set_num_threads(nthreads);
   //#pragma omp for private(j)
   #pragma omp parrallel
   {
      int j, tid = omp_get_thread_num();
      if( tid == 0 ) {
	 nthreads = omp_get_num_threads();
	 printf("%d threads created.\n", nthreads);
      }
      else {
	 #pragma omp critical
	 { i++; }
      }
      fprintf(stdout, "Thread %d started at %f ms.  i = %d\n", tid, (ClockGetTime()-ts)/1.e6, i);
      for(j=0;j<NARR;j++){
	 a[i] += 13.+2.*b[j];
      }
      fprintf(stdout, "Thread %d ended at %f ms.  i = %d\n", tid, (ClockGetTime()-ts)/1.e6, i);
   }
   #pragma omp barrier
   cout<<setprecision(15)<<"All threads done at "<<(ClockGetTime()-ts)/1.e6<<"ms"<<endl;
   delete [] a; a=NULL;
   delete [] b; b=NULL;
   return 0;
}
