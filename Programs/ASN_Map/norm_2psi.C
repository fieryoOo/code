#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
using namespace std;

int calc_azimuth(double lati1, double long1, double lati2, double long2, double *alpha1)
{
  double dlt_lati,dlt_long;
  double Ra,Rb,f;
  double U1,U2;
  double pi;

  pi=4.0*atan(1.0);
  Ra = 6378.137; Rb = 6356.7523142;
  f = 1/298.257223563;
  long1=long1-(int)floor(long1)/360*360;
  if(long1<0) long1+=360;
  long2=long2-(int)floor(long2)/360*360;
  if(long2<0) long2+=360;
  if(lati1==-lati2 && fabs(long1-long2)==180){
  *alpha1=999;
  return 1;
  }
  dlt_lati=fabs(lati2-lati1);
  dlt_long=long2-long1;

if (dlt_long > 180.000)  dlt_long = 360.000000 - dlt_long;
if (dlt_long < -180.000) dlt_long = 360.000 - fabs(dlt_long);
dlt_long = fabs(dlt_long);

U1 = atan((1-f)*tan(lati1/180*pi));
U2 = atan((1-f)*tan(lati2/180*pi));
dlt_long = dlt_long*pi/180;

*alpha1=atan2(cos(U2)*sin(dlt_long), cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(dlt_long))*180/pi;
//*alpha2=atan2(cos(U1)*sin(dlt_long), -sin(U1)*cos(U2) + cos(U1)*sin(U2)*cos(dlt_long))*180/pi;
if( fabs(long2-long1)>180 ) { *alpha1 = 360-*alpha1; }
if( long2 < long1 ) { *alpha1 = 360-*alpha1; }

return 1;
}

#define NSTA 1000
int main (int argc, char *argv[])
{
   FILE *fin, *fout;
   int i,j,k,nsta;
   char buff[300];
   float clong=245.0849, clati=41.415700, longlst[NSTA], latilst[NSTA];
   double inamp[NSTA];
   double A0,A1,A2,phi1,phi2,alpha,pi=3.14159265359;

   if(argc != 7){
      cout<<"Usage: norm_1psi_2psi [input_map] [A0] [A1] [phi1] [A2] [phi2]"<<endl;
      return -1;
   }

   A0=atof(argv[2]);
   A1=atof(argv[3]);
   phi1=atof(argv[4]);
   A2=atof(argv[5]);
   phi2=atof(argv[6]); 

   if((fin=fopen(argv[1],"r"))==NULL) {
      cout<<"Can't open file "<<argv[1]<<endl;
      return 0;
   }
   for(i=0;;i++){
      if((fgets(buff,300,fin))==NULL) break;
      sscanf(buff,"%f %f %lf",&longlst[i],&latilst[i],&inamp[i]);
   }
   nsta=i;
   fclose(fin);

   for(i=0;i<nsta;i++){
      calc_azimuth(clati, clong, latilst[i], longlst[i], &alpha);
      inamp[i] /= A0+A1*sin(pi/180.*(alpha+phi1))+A2*sin(pi/180.*(2*alpha+phi2));
   }

   sprintf(buff,"%s_norm\0",argv[1]);
   fout=fopen(buff,"w");
   for(i=0;i<nsta;i++)
      fprintf(fout,"%f %f %lf\n",longlst[i],latilst[i],inamp[i]);
   fclose(fout);

   return 1;
}

