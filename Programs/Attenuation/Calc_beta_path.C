#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

#define NSTA 2000

double get_dist(double lat1,double lon1,double lat2,double lon2)
{
  double theta,pi,temp;
  double radius=6371;
  pi=4.0*atan(1.0);

  lat1=atan(0.993277*tan(lat1/180*pi))*180/pi;
  lat2=atan(0.993277*tan(lat2/180*pi))*180/pi;

  temp=sin((90-lat1)/180*pi)*cos(lon1/180*pi)*sin((90-lat2)/180*pi)*cos(lon2/180*pi)+sin((90-lat1)/180*pi)*sin(lon1/180*pi)*sin((90-lat2)/180*pi)*sin(lon2/180*pi)+cos((90-lat1)/180*pi)*cos((90-lat2)/180*pi);
  if(temp>1)
    {
      cout<<"warning cos(theta)>1 and correct to 1!!"<<temp<<endl;
      temp=1;
    }
  if(temp<-1)
    {
      cout<<"warning cos(theta)<-1 and correct to -1!!"<<temp<<endl;
      temp=-1;
    }
  theta=fabs(acos(temp));
  return theta*radius;
}

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


struct Site {
   double lon, lat;
   double val;
};

int main (int argc, char *argv[])
{
   if (argc != 7) {
     cout<<"please input [beta_infile] [resolution in km] [source lon] [source lat] [receiver lon] [receiver lat]"<<endl;
     return 0;
   }
 
   FILE *ff;
   char buff[300];
   double rsl = atof(argv[2]), alpha = 0.5/rsl/rsl;
   Site Ss, Sr;
   Ss.lon = atof(argv[3]), Ss.lat = atof(argv[4]);
   Sr.lon = atof(argv[5]), Sr.lat = atof(argv[6]);
   if(Ss.lon<0.) Ss.lon += 360.; if(Sr.lon<0.) Sr.lon += 360.;
  
   if((ff=fopen(argv[1],"r"))==NULL){
     cout<<"Can't open beta file "<<argv[1]<<endl;
     exit(0);
   }
   Site Stmp;
   double diss, disr;
   double sweit = 0., rweit = 0., weight;
   while( fgets(buff, 300, ff)!=NULL ) {
     sscanf(buff,"%lf %lf %lf", &Stmp.lon, &Stmp.lat, &Stmp.val);
     if(Stmp.lon<0) Stmp.lon += 360.;
     //calc_azimuth(clat,clon,lat[i],lon[i],&azi[i]);
     diss = get_dist(Ss.lat, Ss.lon, Stmp.lat, Stmp.lon);
     disr = get_dist(Sr.lat, Sr.lon, Stmp.lat, Stmp.lon);
     if( diss < 3.*rsl ) {
	weight = exp(-alpha*diss*diss);
	Ss.val += weight*Stmp.val;
	sweit += weight;
     }
     if( disr < 3.*rsl ) {
        weight = exp(-alpha*disr*disr);
        Sr.val += weight*Stmp.val;
        rweit += weight;
     }
   }
   fclose(ff);

   if( sweit<1. || rweit<1. ) {
      cout<<"source or receiver location not covered by input beta map"<<endl;
      return 0;
   }

   double dist = get_dist(Ss.lat, Ss.lon, Sr.lat, Sr.lon);
   cout<<Sr.lon<<" "<<Sr.lat<<" "<<log( (Ss.val/sweit) / (Sr.val/rweit) )/dist<<endl;

   return 1;
}
