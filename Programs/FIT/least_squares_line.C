#include<stdio.h>
#include <stdlib.h>
#include<iostream>
#include<math.h>
#include <string>
#include <unistd.h>
#include <float.h>
using namespace std;

#define NBLK 2000

void least_square_fit(int type, double *datx, double *daty, double *sigma, int ndat, double *aout, double *sigmaaout, double *bout, double *sigmabout, double *rms ) {
   if( ndat < 2 ) {
      *aout = *bout = *sigmaaout = *sigmabout = -12345.;
      return;
   }

   int i;
   double a, b, w, x, y, W=0, WX=0, WY=0, WX2=0, WY2=0, WXY=0;
   for(i=0;i<ndat;i++) {
      w = 1./( sigma[i] * sigma[i] );
      x = datx[i]; y = daty[i];
      W += w;
      WX += w * x; WY += w * y;
      WX2 += w * x * x; WY2 += w * y * y;
      WXY += w * x * y; 
   }

   // determine a b accordingly
   if(type == 0) {
      w = 1./(W*WX2-WX*WX);
      a = (W*WXY-WX*WY) * w;
      b = (-WX*WXY+WX2*WY) * w;
   }
   else if(type == 1) {
      w = 1./(W*WXY-WX*WY);
      a=(W*WY2-WY*WY) * w;
      b=(WY*WXY-WY2*WX) * w;
   }
   else {
      cout<<"Line_fit: Wrong input for indep var, stopped!"<<endl;
      exit(0);
   }
   *aout = a; *bout = b;

   if( ndat == 2 ) {
      *sigmaaout = *sigmabout = DBL_MAX;
      return;
   }
   // compute uncertainty in each parameter
   double S2 = 0., k; 
   double dtmp;
   // define k according to ndat from t distribution
   // assuming a 95% conf is equivalent to 2 sigma conf
   if( ndat == 3 ) k = 12.706;
   else if ( ndat == 4 ) k = 4.303;
   else if ( ndat == 5 ) k = 3.182;
   else if ( ndat == 6 ) k = 2.776;
   else k = 1.960+13.8/pow((double)ndat, 1.6); // this could be made better
   k *= 0.5; // now this is 1 sigma
   // compute uncertainties
   for(i=0;i<ndat;i++) {
      dtmp = daty[i] - a * datx[i] - b;
      S2 += dtmp * dtmp;
   }
   S2 /= ndat-2.; *rms = sqrt(S2);
   if( type == 0 ) {
      *sigmaaout = k * sqrt(S2 * W * w);
      *sigmabout = k * sqrt(S2 * WX2 * w);
   }
   else {
      S2 /= a*a;
      *sigmaaout = k * sqrt(S2 * W * w);
      *sigmabout = k * sqrt(S2 * WY2 * w);
   }
}


int main (int argc, char *argv[])
{
   if(argc != 3){
      printf("Usage: least_squares_line.C [input file (x y sigma)] [indep var (0 for x, 1 for y)]\n");
      exit(-1);
   }

   FILE *ff;
   char buff[300];
   int i, itmp, ndat;
   //double datx[NDAT], daty[NDAT], sigma[NDAT];
	double *datx = NULL, *daty = NULL, *sigma = NULL;

   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Can't open file: "<<argv[1]<<endl;
      exit(0);
   }

	int nsize = 0;
   for(i=0; fgets(buff, 300, ff); ) {
		if( i >= nsize ) {
			nsize += NBLK;
			datx = (double*) realloc (datx, nsize * sizeof(double));
			daty = (double*) realloc (daty, nsize * sizeof(double));
			sigma = (double*) realloc (sigma, nsize * sizeof(double));
		}
      itmp = sscanf(buff,"%lf %lf %lf", &datx[i], &daty[i], &sigma[i]);
      if( itmp == 2 ) sigma[i]=1.;
      else if (itmp != 3) continue;
      i++;
   }
   fclose(ff);
   ndat=i;

   double a, b, sigmaa, sigmab, rms;
   least_square_fit(atoi(argv[2]), datx, daty, sigma, ndat, &a, &sigmaa, &b, &sigmab, &rms);
   cout<<a<<" "<<b<<" "<<sigmaa<<" "<<sigmab<<" "<<rms<<std::endl;

	free(datx); free(daty); free(sigma);

   return 1;
}
