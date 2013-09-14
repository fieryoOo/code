#include<stdio.h>
#include <stdlib.h>
#include<iostream>
#include<math.h>
#include <string>
#include <unistd.h>
using namespace std;

#define NDAT 2000

int main (int argc, char *argv[])
{
   if(argc != 2){
      printf("Usage: least_squares_parabola.C [input file (x y weight)]\n");
      exit(-1);
   }

   FILE *ff;
   char buff[300];
   int i, ndat;
   double datx[NDAT], daty[NDAT], weit[NDAT], a, b, c, std, stdw;

   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Can't open file: "<<argv[1]<<endl;
      exit(0);
   }
   for(i=0;;i++) {
      if(fgets(buff,300,ff)==NULL) break;
      sscanf(buff,"%lf %lf %lf",&datx[i],&daty[i],&weit[i]);
      weit[i] = 1./weit[i];
   }
   fclose(ff);
   ndat=i;

   double X=0, Y=0, X2=0, Y2=0, XY=0;
   for(i=0;i<ndat;i++) {
      X += datx[i];
      Y += daty[i];
      X2 += pow(datx[i],2);
      Y2 += pow(daty[i],2);
      XY += datx[i]*daty[i];
   }

   double varX, varY, cov;
   varX = (X2-X*X/ndat)/(ndat-3);
   varY = (Y2-Y*Y/ndat)/(ndat-3);
   cov = (XY-X*Y/ndat)/(ndat-3);

   cout<<"Covariance matrix:  "<<endl;
   cout<<varX<<"  "<<cov<<endl;
   cout<<cov<<"  "<<varY<<endl;

   return 1;
}
