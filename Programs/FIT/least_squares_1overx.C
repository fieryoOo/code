#include<stdio.h>
#include <stdlib.h>
#include<iostream>
#include<math.h>
#include <string>
#include <unistd.h>
using namespace std;

#define NDAT 100000

int main (int argc, char *argv[])
{
   if(argc != 3){
      printf("Usage: least_squares_line.C [input file] [indep var (0 for x, 1 for y)]\n");
      exit(-1);
   }

   FILE *ff;
   char buff[300];
   int i, ndat;
   double datx[NDAT], daty[NDAT], std;

   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Can't open file: "<<argv[1]<<endl;
      exit(0);
   }
   for(i=0;;i++) {
      if(fgets(buff,300,ff)==NULL) break;
      sscanf(buff,"%lf %lf",&datx[i],&daty[i]);
      datx[i] = 1./datx[i];
   }
   fclose(ff);
   ndat=i;

   double a, b, X=0, Y=0, X2=0, Y2=0, XY=0;
   for(i=0;i<ndat;i++) {
      X += datx[i];
      Y += daty[i];
      X2 += pow(datx[i],2);
      Y2 += pow(daty[i],2);
      XY += datx[i]*daty[i];
   }
   std=0;
   for(i=0;i<ndat;i++) std += pow(daty[i]-a*datx[i]-b,2);
   if(atoi(argv[2])==0) {
      a = (ndat*XY-X*Y)/(ndat*X2-X*X);
      b = (-X*XY+X2*Y)/(ndat*X2-X*X);
   }
   else if(atoi(argv[2])==1) {
      a=(ndat*Y2-Y*Y)/(ndat*XY-X*Y);
      b=(Y*XY-Y2*X)/(ndat*XY-X*Y);
      std /= a*a;
   }
   else {
      cout<<"Line_fit: Wrong input for indep var, stopped!"<<endl;
      exit(0);
   }
   std = sqrt(std/(ndat-1));

   cout<<a<<" "<<b<<" "<<std<<endl;

   return 1;
}
