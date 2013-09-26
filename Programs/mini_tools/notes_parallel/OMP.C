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
   int i=1000;
   float *a = new float[NARR], *b = new float[NARR];
   double ts = ClockGetTime();
   //omp_set_num_threads(20);
   int nthreads = omp_get_num_threads();
   //#pragma omp parallel for private(j)
   #pragma omp parallel shared(i)
   {
      int j;
      printf("Thread %d start at %lf ms.  i = %d\n", omp_get_thread_num(), (ClockGetTime()-ts)/1.e6, i);
      for(j=0;j<NARR;j++){
	 a[i] += 13.+2.*b[j];
      }
      #pragma omp critical
      i++;
      printf("Thread %d end at %lf ms.  i = %d\n", omp_get_thread_num(), (ClockGetTime()-ts)/1.e6, i);
   }
   #pragma omp barrier
   cout<<setprecision(15)<<"All "<<nthreads<<" threads done at "<<(ClockGetTime()-ts)/1.e6<<"ms"<<endl;
   return 0;
}
