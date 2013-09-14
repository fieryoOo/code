#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
using namespace std;

#define NDAT 10000

int main(int argc, char *argv[])
{
   if(argc!=2){
      cout<<"Usage: Correlation [in_file (col1: dat1  col2: dat2)]"<<endl;
      return 0;
   }

   FILE *ff;
   double dat1[NDAT], dat2[NDAT], avg1, avg2, corr, weight1, weight2;
   int i, ndat;
   char buff[300];

   avg1=0.; avg2=0.;
   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Can't open file "<<argv[1]<<" to read"<<endl;
      return 0;
   }
   for(i=0;;i++) {
      if((fgets(buff, 300, ff))==NULL) break;
      sscanf(buff, "%lf %lf", &dat1[i], &dat2[i]);
      avg1 += dat1[i];
      avg2 += dat2[i];
   }
   fclose(ff);
   ndat = i;
   avg1 /= ndat;
   avg2 /= ndat;

   weight1 = 0.; weight2 = 0.;
   for(i=0;i<ndat;i++) {
      corr+=(dat1[i]-avg1)*(dat2[i]-avg2);
      weight1+=pow(dat1[i]-avg1,2);
      weight2+=pow(dat2[i]-avg2,2);
   }
   corr=corr/sqrt(weight1*weight2);
   cout<<corr<<endl;

   return 1;
}
