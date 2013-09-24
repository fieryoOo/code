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
   int i, j;
   float *a = new float[NARR], *b = new float[NARR];
   double ts = ClockGetTime();
   omp_set_num_threads(10);
   #pragma omp for private j
   for(i=0;i<10;i++) {
      cout<<setprecision(15)<<"Thread "<<i<<"start at "<<(ClockGetTime()-ts)/1.e6<<"ms"<<endl;
      for(j=0;j<NARR;j++){
	 a[i] += 13.+2.*b[j];
      }
   }
   cout<<setprecision(15)<<"All threads done at "<<(ClockGetTime()-ts)/1.e6<<"ms"<<endl;
   return 0;
}
