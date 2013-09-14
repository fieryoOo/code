#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <string>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
using namespace std;

//--------Clock, swap & median----------//
uint64_t ClockGetTime()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
}


#define SBUF 10000

int main (int argc, char *argv[])
{
   if(argc != 2){
      printf("Usage: %s [#operations]\n", argv[0]);
      exit(-1);
   }

   // variables
   double a = 3.14, b = 108.96, c;
   //warm up
   int i, nop = atoi(argv[1]);
   for(i=0;i<10000;i++) {
      c = a+b;
      c = a-b;
      c = a*b;
      c = a/b;
      c = sin(b);
   }

   double ts, te;
   // summation
   ts=ClockGetTime();
   for(i=0; i<nop; i++) c = a+b;
   te=ClockGetTime();
   cout<<nop<<" summations:\t"<<te-ts<<"msec "<<endl;
   // substraction
   ts=ClockGetTime();
   for(i=0; i<nop; i++) c = a-b;
   te=ClockGetTime();
   cout<<nop<<" substraction:\t"<<te-ts<<"msec "<<endl;
   // multiplication
   ts=ClockGetTime();
   for(i=0; i<nop; i++) c = a*b;
   te=ClockGetTime();
   cout<<nop<<" multiplications:\t"<<te-ts<<"msec "<<endl;
   // divide by 2
   ts=ClockGetTime();
   for(i=0; i<nop; i++) c = a/2;
   te=ClockGetTime();
   cout<<nop<<" divisions(by 2):\t"<<te-ts<<"msec "<<endl;
   // division
   ts=ClockGetTime();
   for(i=0; i<nop; i++) c = a/b;
   te=ClockGetTime();
   cout<<nop<<" divisions:\t"<<te-ts<<"msec "<<endl;
   // sin(x)
   ts=ClockGetTime();
   for(i=0; i<nop; i++) c = sin(b);
   te=ClockGetTime();
   cout<<nop<<" sin(x)s:\t\t"<<te-ts<<"msec "<<endl;

   return 1;
}
