#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

int discard_by_std (double *dat, int ndat, double std_coef, int *count)
{
   int i, ndat2;
   double avg, std;

   avg=0;
   for(i=0;i<ndat;i++) avg += dat[i];
   avg /= ndat; 
   std=0;
   for(i=0;i<ndat;i++) std += pow(dat[i]-avg,2);
   std = sqrt(std/(ndat-1));
   for(i=0;i<ndat;i++) if(fabs(dat[i]-avg)>std_coef*std) dat[i]=-1;

   avg=0; ndat2=0;
   for(i=0;i<ndat;i++) {
      if(dat[i]==-1) continue;
      avg += dat[i];
      ndat2++;
   }
   avg /= ndat2;
   std=0;
   for(i=0;i<ndat;i++) {
      if(dat[i]==-1) continue;
      std += pow(dat[i]-avg,2);
   }
   std = sqrt(std/(ndat2-1));
   for(i=0;i<ndat;i++) if(fabs(dat[i]-avg)>std_coef*std) dat[i]=-1;

   *count = 0;
   for(i=0;i<ndat;i++) if(dat[i]!=-1) *count = *count + 1;

   return 1;
}
