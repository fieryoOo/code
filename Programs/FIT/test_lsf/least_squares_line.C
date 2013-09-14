#include<stdio.h>
#include <stdlib.h>
#include<iostream>
#include<math.h>
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
   if(argc != 3){
      printf("Usage: least_squares_line.C [input file] [indep var (0 for x, 1 for y)]\n");
      exit(-1);
   }

   FILE *ff;
   char buff[300];
   int i, ndat, nbuf;
   double *datx=NULL, *daty=NULL, std;

   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Can't open file: "<<argv[1]<<endl;
      exit(0);
   }
   for(nbuf=0,i=0;;i++) {
      if(fgets(buff,300,ff)==NULL) break;
      if(nbuf*SBUF-1<i) {
	 nbuf++;
	 datx = (double *)realloc(datx, nbuf*SBUF*sizeof(double));
	 daty = (double *)realloc(daty, nbuf*SBUF*sizeof(double));
      }
      sscanf(buff,"%lf %lf",&datx[i],&daty[i]);
   }
   fclose(ff);
   ndat=i;

   double ts=ClockGetTime();
   double a, b, X=0, Y=0, X2=0, Y2=0, XY=0;
   for(i=0;i<ndat;i++) {
      X += datx[i];
      Y += daty[i];
      X2 += pow(datx[i],2);
      Y2 += pow(daty[i],2);
      XY += datx[i]*daty[i];
   }
   if(atoi(argv[2])==0) {
      a = (ndat*XY-X*Y)/(ndat*X2-X*X);
      b = (-X*XY+X2*Y)/(ndat*X2-X*X);
   }
   else if(atoi(argv[2])==1) {
      a=(ndat*Y2-Y*Y)/(ndat*XY-X*Y);
      b=(Y*XY-Y2*X)/(ndat*XY-X*Y);
   }
   else {
      cout<<"Line_fit: Wrong input for indep var, stopped!"<<endl;
      exit(0);
   }
   double te=ClockGetTime();
   std=0;
   for(i=0;i<ndat;i++) std += pow(daty[i]-a*datx[i]-b,2);
   if(atoi(argv[2])==1) std /= a*a;
   std = sqrt(std/(ndat-1));
   free(datx); free(daty);

   cout<<a<<" "<<b<<" "<<te-ts<<"msec "<<endl;

   return 1;
}
