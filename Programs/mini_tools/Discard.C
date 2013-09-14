#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

#define NMOD 5000

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
   cout<<"original:  mean "<<avg<<"  std "<<std<<endl;
   for(i=0;i<ndat;i++) if(fabs(dat[i]-avg)>std_coef*1.5*std) dat[i]=-1;

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
   cout<<"After 1st round:  mean "<<avg<<"  std "<<std<<endl;
   for(i=0;i<ndat;i++) if(fabs(dat[i]-avg)>std_coef*std) dat[i]=-1;

   *count = 0; avg=0;
   for(i=0;i<ndat;i++) {
       if(dat[i]==-1) continue;
       *count = *count + 1;
       avg += dat[i];
   }
   avg /= *count; std=0;
   for(i=0;i<ndat;i++) {
      if(dat[i]==-1) continue;
      std += pow(dat[i]-avg,2);
   }
   std = sqrt(std/(*count-1));
   cout<<"After 2nd round:  mean "<<avg<<"  std "<<std<<endl;

   return 1;
}

int main (int argc, char *argv[])
{
   if(argc != 3) {
      cout<<"Usage: discard [in_file] [std_coef]"<<endl;
      exit(-1);
   }

   FILE *ff;
   int i, ndat, nmod = 0;
   char buff[300];
   double *datx = NULL, *daty = NULL;

   if((ff=fopen(argv[1], "r")) == NULL) {
      cout<<"Can't open file: "<<argv[1]<<endl;
      exit(0);
   }
   for(i=0;;i++) {
      if(nmod*NMOD<=i) {
         datx = (double *) realloc( datx, (++nmod)*NMOD * sizeof(double));
	 daty = (double *) realloc( daty, nmod*NMOD * sizeof(double));
      }
      if(fgets(buff, 300, ff) == NULL) break;
      sscanf(buff, "%lf %lf", &datx[i], &daty[i]);
   }
   fclose(ff);
   ndat=i;

   int count;
   discard_by_std (&daty[0], ndat, atof(argv[2]), &count);

   sprintf(buff, "%s_dscd", argv[1]);
   ff=fopen(buff, "w");
   for(i=0;i<ndat;i++) {
      if(daty[i]==-1) continue;
      fprintf(ff, "%lf %lf\n", datx[i], daty[i]);
   }
   fclose(ff);
   free(datx); free(daty);

   return 1;
}
