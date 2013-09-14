#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <time.h>
using namespace std;

#define DGR 0.1

int calc_azimuth(double lati1, double long1, double lati2, double long2, double *alpha2)
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
  *alpha2=999;
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

//*alpha1=atan2(cos(U2)*sin(dlt_long), cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(dlt_long))*180/pi;
*alpha2=atan2(cos(U1)*sin(dlt_long), -sin(U1)*cos(U2) + cos(U1)*sin(U2)*cos(dlt_long))*180/pi;
if( fabs(long2-long1)>180 ) { *alpha2 = 360-*alpha2; }
if( long2 < long1 ) { *alpha2 = 360-*alpha2; }

return 1;
}

int calc_dist(double lati1, double long1, double lati2, double long2, double *dist)
{
  int i,flag=0;
  double latio1,longo1,dlt_lati,dlt_long;
  double Ra,Rb,R,Rx,f;
  double U1,U2;
  double Ds,alpha2;
  double ctr_angl;
  double pi;
  double cv,cv1,cv2,cv3,cv4,cv5,cvC,numda1;
  double numda;
  double mius, cvA, cvB, deltacv;

  pi=4.0*atan(1.0);
  Ra = 6378.137; Rb = 6356.7523142;
  f = 1/298.257223563;
begin:
  long1=long1-(int)floor(long1)/360*360;
  if(long1<0) long1+=360;
  long2=long2-(int)floor(long2)/360*360;
  if(long2<0) long2+=360;
  if(lati1==-lati2 && fabs(long1-long2)==180){
  *dist=20003.917357;
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

numda = dlt_long;
numda1 = numda;
i=0;
do {
i++;
numda = numda1;
cv1 =  sqrt( (cos(U2)*sin(numda))*(cos(U2)*sin(numda))+ (cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda))*(cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda)) );
cv2 = sin(U1)*sin(U2)+ cos(U1)*cos(U2)*cos(numda);
cv = atan2(cv1,cv2);
if(cv==0)cv+=0.00000000001;
cv3 = cos(U1)*cos(U2)*sin(numda)/sin(cv);
cv4 = 1 - cv3*cv3;
if(cv4==0)cv4+=0.00000000001;
cv5 = cos(cv) - 2*sin(U1)*sin(U2)/cv4;
cvC = f/16*cv4*(4 + f*(4 - 3*cv4));
numda1 = dlt_long + (1-cvC)*f*cv3*(cv + cvC*cv1*(cv5 + cvC*cv2*(-1 +2*cv5*cv5)));
if(i>50){
flag=1;
latio1=lati1*pi/180.;
Rx=sqrt((pow(Ra*Ra*cos(latio1),2)+pow(Rb*Rb*sin(latio1),2))/(pow(Ra*cos(latio1),2)+pow(Rb*sin(latio1),2)));
Rx=pi*(Rx+Ra)/2.;
latio1=lati1;
longo1=long1;
lati1=lati2;
long1=long2;
lati2=-latio1;
long2=longo1+180;
goto begin;
}
} while (fabs(numda - numda1) > 1e-10);
mius = cv4*(Ra*Ra - Rb*Rb)/(Rb*Rb);
cvA = 1+mius/16384*(4096 + mius*(-768 + mius*(320 - 175*mius)));
cvB = mius/1024*(256+ mius*(-128 + mius*(74 - 47*mius)));
deltacv = cvB*cv1*(cv5 +cvB/4*(cv2*(-1 + 2*cv5*cv5)-cvB/6*cv5*(-3+4*cv1*cv1)*(-3+4*cv5*cv5) ));
*dist = Rb * cvA *(cv - deltacv);

if(flag==1){
alpha2=atan2(cos(U1)*sin(dlt_long), -sin(U1)*cos(U2) + cos(U1)*sin(U2)*cos(dlt_long))*180/pi;
if( fabs(long2-long1)>180 ) { alpha2 = 360-alpha2;}
if( long2 < long1 ) { alpha2 = 360-alpha2; }
if(alpha2>180) alpha2=360-alpha2;
alpha2=fabs(90-alpha2);
Ds=Rx*(90-alpha2)/90.+(Ra+Rb)*pi*alpha2/180.;
*dist=Ds-*dist;
}

return 1;
}

int   calc_grd_gcp(double clong, double clati, char *file, double longmin, int npts_long, double latimin, int npts_lati, double **dat, double **gradr)
{
   FILE *ff;
   char buff[300], buff2[300];
   double longtmp, latitmp;
   double radius=6371.1391285, pi=4.0*atan(1.0), dgr=0.05, dattmp;
   double grdtx, grdty, dis, alpha, weight;
   int i, j, ii, jj, ista, nrand, ntmpp;
   double dislong[npts_lati], dislati, distmp;
   double gradtr, grdtmp, gradx, grady;

   for(i=0;i<npts_long;i++) for(j=0;j<npts_lati;j++){
      gradr[i][j] = 1e10;
      dat[i][j] = -1;
   }
   for(i=1;i<npts_lati-1;i++){
      dislati=atan(0.993277*tan((latimin+i*dgr)/180*pi))*180/pi;
      dislong[i]=radius*sin((90-dislati)/180*pi)*dgr/180*pi;
   }
   dislati=radius*dgr/180*pi;

   if((ff=fopen(file,"r"))==NULL){
     printf("Can't open %s to read!\n",file);
     return 0;
   }
   for(;;){
      if(fgets(buff, 300, ff) == NULL ) break;
      if((sscanf(buff,"%lf %lf %lf", &longtmp, &latitmp, &dattmp))!=3) {
         cout<<"Wrong format! Stopped: "<<buff<<endl;
         break;
      }
      i=int((longtmp-longmin)/DGR+0.1);
      j=int((latitmp-latimin)/DGR+0.1);
      if(i<0||i>=npts_long||j<0||j>=npts_lati) continue;
      dat[i][j]=dattmp;
   }
   fclose(ff);

   for(i=1;i<npts_long-1;i++)
      for(j=1;j<npts_lati-1;j++){
         if( dat[i-1][j-1] == -1 || dat[i][j-1] == -1 || dat[i+1][j-1] == -1 || dat[i-1][j] == -1 || dat[i+1][j] == -1 || dat[i-1][j+1] == -1 || dat[i][j+1] == -1 || dat[i+1][j+1] == -1 ) continue;
         gradx=(dat[i+1][j]-dat[i-1][j])/2.0/dislong[j];
         grady=(dat[i][j+1]-dat[i][j-1])/2.0/dislati;
         gradtr=sqrt(gradx*gradx+grady*grady);
         distmp=sqrt(dislong[j]*dislong[j]+dislati*dislati);
         grdtx=(dat[i+1][j-1]-dat[i-1][j+1])/2.0/distmp;
         grdty=(dat[i+1][j+1]-dat[i-1][j-1])/2.0/distmp;
         alpha=2*atan(dislati/dislong[j]);
         alpha=atan((grdtx/grdty-cos(alpha))/sin(alpha));
         grdtmp=fabs(grdty/cos(alpha));
         //cout<<gradtr[i][j]<<" "<<grdtmp<<endl;
         if(fabs((gradtr-grdtmp)/gradtr)>0.15) gradr[i][j]=1e10;
         else {
           longtmp=longmin+i*DGR; latitmp=latimin+j*DGR;
           alpha=0;
           calc_azimuth(clati, clong, latitmp, longtmp, &alpha);
           alpha=alpha*pi/180.;
           gradr[i][j]=gradx*sin(alpha)+grady*cos(alpha);
         }
      }


   return 1;
}


int main (int argc, char *argv[])
{
   if(argc!=3){
      cout<<"usage:csta_map_gradient [input_beta_map] [input_pvel_map]"<<endl;
      return 0;
   }
   FILE *ff1, *ff2;
   char fname[100], *tmpc;
   double blong[2000], blati[2000], dist[2000];
   double clon=245.124300, clat=41.100300;
   double lonmin=235, latmin=25, lon, lat;
   int i,j;
   int nptslon=int((265-lonmin)/DGR)+1, nptslat=int((50-latmin)/DGR)+1;
   double *beta[nptslon], *grdb[nptslon], *pvel[nptslon], *grdc[nptslon];
//   double *beta[nptslon], *grdxb[nptslon], *grdyb[nptslon], *pvel[nptslon], *grdxc[nptslon], *grdyc[nptslon];

   for(i=0;i<nptslon;i++){
      beta[i] = (double *) malloc( nptslat * sizeof(double));
      grdb[i] = (double *) malloc( nptslat * sizeof(double));
//      grdxb[i] = (double *) malloc( nptslat * sizeof(double));
//      grdyb[i] = (double *) malloc( nptslat * sizeof(double));
      pvel[i] = (double *) malloc( nptslat * sizeof(double));
      grdc[i] = (double *) malloc( nptslat * sizeof(double));
//      grdxc[i] = (double *) malloc( nptslat * sizeof(double));
//      grdyc[i] = (double *) malloc( nptslat * sizeof(double));      
   }
//   double beta[nptslon][nptslat], grdxb[nptslon][nptslat], grdyb[nptslon][nptslat], pvel[nptslon][nptslat], grdxc[nptslon][nptslat], grdyc[nptslon][nptslat];

   calc_grd_gcp(clon, clat, argv[1], lonmin, nptslon, latmin, nptslat, beta, grdb);
//for(i=0;i<nptslon;i++)for(j=0;j<nptslat;j++) if(beta[i][j]!=-1)
  // cout<<beta[i][j]<<" "<<grdxb[i][j]<<" "<<grdyb[i][j]<<endl;
   calc_grd_gcp(clon, clat, argv[2], lonmin, nptslon, latmin, nptslat, pvel, grdc);
//for(i=0;i<nptslon;i++)for(j=0;j<nptslat;j++) if(pvel[i][j]!=-1)
  // cout<<pvel[i][j]<<" "<<grdxc[i][j]<<" "<<grdyc[i][j]<<endl;

   double dec[nptslon][nptslat];
   double site[nptslon][nptslat];
   double prp[nptslon][nptslat];
   for(i=0;i<nptslon;i++){
      for(j=0;j<nptslat;j++){
         if(grdb[i][j]>=0.9e10||grdc[i][j]>=0.9e10) continue;
         site[i][j]=-grdb[i][j]/beta[i][j];
         //sitey[i][j]=-grdyb[i][j]/beta[i][j];
         prp[i][j]=-0.5*grdc[i][j]/pvel[i][j];
         //prpy[i][j]=-0.5*grdyc[i][j]/pvel[i][j];
         dec[i][j]=site[i][j]+prp[i][j];
         //decy[i][j]=sitey[i][j]+prpy[i][j];
      }
   }

   sprintf(fname,"alpha_beta_g_map",tmpc);
   ff1=fopen(fname,"w");
//   sprintf(fname,"alpha_betay_map",tmpc);
//   ff2=fopen(fname,"w");
   for(i=0;i<nptslon;i++)
      for(j=0;j<nptslat;j++){
         if(grdb[i][j]>=0.9e10||grdc[i][j]>=0.9e10) continue;
         lon=lonmin+i*DGR; lat=latmin+j*DGR;
         fprintf(ff1,"%8.4f  %8.4f  %8g\n", lon, lat, site[i][j]);
//         fprintf(ff2,"%8.4f  %8.4f  %8g\n", lon, lat, site[i][j]);
      }
   fclose(ff1);
//   fclose(ff2);

   sprintf(fname,"alpha_pvel_g_map",tmpc);
   ff1=fopen(fname,"w");
//   sprintf(fname,"alpha_pvely_map",tmpc);
//   ff2=fopen(fname,"w");
   for(i=0;i<nptslon;i++)
      for(j=0;j<nptslat;j++){
         if(grdb[i][j]>=0.9e10||grdc[i][j]>=0.9e10) continue;
         lon=lonmin+i*DGR; lat=latmin+j*DGR;
         fprintf(ff1,"%8.4f  %8.4f  %8g\n", lon, lat, prp[i][j]);
//         fprintf(ff2,"%8.4f  %8.4f  %8g\n", lon, lat, prp[i][j]);
      }
   fclose(ff1);
//   fclose(ff2);

   sprintf(fname,"alpha_g_map",tmpc);
   ff1=fopen(fname,"w");
//   sprintf(fname,"alpha_y_map",tmpc);
//   ff2=fopen(fname,"w");
   for(i=0;i<nptslon;i++)
      for(j=0;j<nptslat;j++){
         if(grdb[i][j]>=0.9e10||grdc[i][j]>=0.9e10) continue;
         lon=lonmin+i*DGR; lat=latmin+j*DGR;
         fprintf(ff1,"%8.4f  %8.4f  %8g\n", lon, lat, dec[i][j]);
//         fprintf(ff2,"%8.4f  %8.4f  %8g\n", lon, lat, dec[i][j]);
      }
   fclose(ff1);
//   fclose(ff2);


   return 1;
}
