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


int main (int argc, char *argv[])
{
  if (argc != 4) {
    cout<<"please input [alpha_file] [clon] [clat]"<<endl;
    return 0;
  }
 
  FILE *ff;
  char buff[300];
  int i, j, nsta;
  double dat[NSTA], lon[NSTA], lat[NSTA], azi[NSTA], dis[NSTA], bin_num[NSTA];
  double clon=atof(argv[2]), clat=atof(argv[3]);
  if(clon<0) clon += 360.;
  
  if((ff=fopen(argv[1],"r"))==NULL){
    cout<<"Can't open alpha file "<<argv[1]<<endl;
  } 
  for(i=0;;i++) {
    if(fgets(buff, 300, ff)==NULL) break;
    sscanf(buff,"%lf %lf %lf", &lon[i], &lat[i], &dat[i]);
    if(lon[i]<0) lon[i] += 360.;
    calc_azimuth(clat,clon,lat[i],lon[i],&azi[i]);
    dis[i]=get_dist(clat, clon, lat[i], lon[i]);
  }
  fclose(ff);
  nsta=i;
  int bin=20, nbin=18, nbine=0, azicv, nesta, ndat0[nbin], binsta=20, ibin, iazi, idat;
   for(ibin=0;ibin<nbin;ibin++) ndat0[ibin]=0;
   for(i=0;i<nsta;i++) {
      ibin=(int)floor((azi[i]+bin/2.)/bin);
      if(ibin>=nbin) ibin=0;
      ndat0[ibin]++;
   }
   nesta=0;
   for(ibin=0;ibin<nbin;ibin++) if(ndat0[ibin]>binsta/3.) { nbine++; nesta+=ndat0[ibin]; }
   azicv=(nbine-1)*bin;
   int bin2=4, nbin2=90, ndat1[nbin2];
   if(nbine<3) {
      for(ibin=0;ibin<nbin2;ibin++) ndat1[ibin]=0;
      for(i=0;i<nsta;i++) {
         ibin=(int)floor((azi[i]+bin2/2.)/bin2);
         if(ibin>=nbin2) ibin=0;
         ndat1[ibin]++;
      }
      nbine=0; nesta=0;
      for(ibin=0;ibin<nbin;ibin++) {
         if(ndat0[ibin]<=binsta/3.) continue;
         for(i=ibin*bin/bin2-2;i<=ibin*bin/bin2+2;i++) {
            j=i;
            if(j<0) j+=nbin2; else if(j>=nbin2) j-=nbin2;
            if(ndat1[j]>=binsta/5.) { nbine++; nesta+=ndat1[j]; }
         }
      }
      azicv=(nbine-1)*bin2;
   }
   bin=(int)floor(azicv/(float(nesta)/binsta)/2.+0.5);
   if(bin==1) bin=2;
   if(bin>20) { cout<<": Stations not dense enough. Skipped!"<<endl; return 0; }
   if(bin == 7 || bin == 9 || bin == 11 || bin == 13 || bin == 16 || bin == 19) bin -= 1;
   if(bin == 14 || bin == 17) bin += 1;
   nbin=(int)(359.9/bin)+1;
   cout<<"azicv: "<<azicv<<"  nesta: "<<nesta<<"  bin: "<<bin<<"  nbin: "<<nbin<<endl;

   sprintf(buff,"mkdir -p interference_%s", argv[1]);
   system(buff);
   double tazi, datx[500], daty[500];
   for(ibin=0;ibin<nbin;ibin++){
      iazi=ibin*bin;
      idat=0;
      for(i=0;i<nsta;i++) {
         tazi=azi[i];
         if( (tazi >= iazi+bin && tazi < iazi-bin+360.) || (tazi < iazi-bin && tazi >= iazi+bin-360.)) continue;
         datx[idat] = dis[i];
         daty[idat] = dat[i]*sqrt(dis[i]);
         idat++;
      }
      if(idat<10) continue;
      sprintf(buff,"interference_%s/%s_%dazi_%dhbin",argv[1],argv[1],iazi,bin);
      ff=fopen(buff,"w");
      for(i=0;i<idat;i++) fprintf(ff,"%lf %lf\n", datx[i], daty[i]);
      fclose(ff);
   }

  return 1;
}
