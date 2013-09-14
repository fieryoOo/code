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


void hex2bin( char* hex, short* bin )
{
   short i, j, hv, len = strlen(hex);

   for(i=0;i<len;i++){
      if (hex[i] >= '0' && hex[i] <= '9') hv = hex[i] - '0';
      else hv = hex[i] - 'A' + 10;
      for(j=i*4+3;j>i*4;j--){
         bin[j]= hv % 2;
         hv = hv/2;
      }
      bin[i*4]=hv;
   }
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

int calc_grdt(char *file, double *slong, double *slati, double *sdat, double *grdt, double *azi, int *nsta)
{
   FILE *ff;
   char buff[300], buff2[300];
   float longmin=360, longmax=0, latimin=90, latimax=-90, longtmp, latitmp;
   double radius=6371.1391285, pi=4.0*atan(1.0), dgr=0.05, dattmp;
   double grdtx, grdty, dis, alpha, weight;
   int i, j, ii, jj, ista, nrand, ntmpp;

   if((ff=fopen(file,"r"))==NULL){
     printf("Can't open %s to read!\n",file);
     return 0;
   }
   for(*nsta=0;;){
      if(fgets(buff, 300, ff) == NULL ) break;
      if((sscanf(buff,"%lf %lf %lf", &slong[*nsta], &slati[*nsta], &sdat[*nsta]))!=3) {
         cout<<"Wrong format! Stopped: "<<buff<<endl;
         break;
      }
      if(longmin>slong[*nsta]) longmin=slong[*nsta];
      if(longmax<slong[*nsta]) longmax=slong[*nsta];
      if(latimin>slati[*nsta]) latimin=slati[*nsta];
      if(latimax<slati[*nsta]) latimax=slati[*nsta];
      *nsta=*nsta+1;
   }
   fclose(ff);
   longmin=floor(longmin/dgr)*dgr;
   latimin=floor(latimin/dgr)*dgr;
   longmax=ceil(longmax/dgr)*dgr;
   latimax=ceil(latimax/dgr)*dgr;
   int npts_long=int((longmax-longmin)/dgr+1);
   int npts_lati=int((latimax-latimin)/dgr+1);
   double dislong[npts_lati], dislati, distmp;
   double gradtr[npts_long][npts_lati], grdtmp, azig[npts_long][npts_lati], dat[npts_long][npts_lati];
   for(i=0;i<npts_long;i++) for(j=0;j<npts_lati;j++){
      gradtr[i][j] = 1e10;
      dat[i][j] = -1;
   }
   for(i=1;i<npts_lati-1;i++){
      dislati=atan(0.993277*tan((latimin+i*dgr)/180*pi))*180/pi;
      dislong[i]=radius*sin((90-dislati)/180*pi)*dgr/180*pi;
   }
   dislati=radius*dgr/180*pi;

   nrand=rand();
   sprintf(buff,"region_%d",nrand);
   ff=fopen(buff,"w");
   fprintf(ff, "-R%d/%d/%d/%d\0", (int)floor(longmin), (int)ceil(longmax), (int)floor(latimin), (int)ceil(latimax));
   fclose(ff);
   sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel_positive %s region_%d %f", file, nrand, dgr);
   system(buff);

   sprintf(buff,"%s.HD",file);
   if((ff=fopen(buff,"r"))==NULL){
      printf("Can't open %s to read!\n",buff);
      return 0;
   }
   for(;;){
      if( fgets(buff, 300, ff) == NULL ) break;
      sscanf(buff,"%f %f %lf", &longtmp, &latitmp, &dattmp);
      i=int((longtmp-longmin)/dgr+0.1);
      j=int((latitmp-latimin)/dgr+0.1);
      if(i<0||i>=npts_long||j<0||j>=npts_lati) continue;
      dat[i][j]=dattmp;
   }
   fclose(ff);

   for(ista=0;ista<*nsta;ista++){
      ii=int((slong[ista]-longmin)/dgr);
      jj=int((slati[ista]-latimin)/dgr);
      for(i=ii;i<ii+2;i++)
         for(j=jj;j<jj+2;j++){
            if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
            if(gradtr[i][j] != 1e10) continue;
            grdtx=(dat[i+1][j]-dat[i-1][j])/2.0/dislong[j];
            grdty=(dat[i][j+1]-dat[i][j-1])/2.0/dislati;
            gradtr[i][j]=sqrt(grdtx*grdtx+grdty*grdty);
            azig[i][j]=pi/2-atan2(grdty,grdtx);
            if(azig[i][j]<0-1e-8)azig[i][j]+=2*pi;
            distmp=sqrt(dislong[j]*dislong[j]+dislati*dislati);
            grdtx=(dat[i+1][j-1]-dat[i-1][j+1])/2.0/distmp;
            grdty=(dat[i+1][j+1]-dat[i-1][j-1])/2.0/distmp;
            alpha=2*atan(dislati/dislong[j]);
            alpha=atan((grdtx/grdty-cos(alpha))/sin(alpha));
            grdtmp=fabs(grdty/cos(alpha));
            //cout<<gradtr[i][j]<<" "<<grdtmp<<endl;
            if(fabs((gradtr[i][j]-grdtmp)/gradtr[i][j])>0.15) gradtr[i][j]=1e20;
            else gradtr[i][j]=(gradtr[i][j]+grdtmp)/2.;
         }
   }

   for(ista=0;ista<*nsta;ista++){
      ii=int((slong[ista]-longmin)/dgr);
      jj=int((slati[ista]-latimin)/dgr);
      ntmpp=0; weight=0; grdt[ista]=0; azi[ista]=0;
      for(i=ii;i<ii+2;i++)
         for(j=jj;j<jj+2;j++){
            if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
            if(gradtr[i][j]>=1e10) continue;
            ntmpp++;
            dis=sqrt(pow((longmin+i*dgr-slong[ista])*dislong[ista],2)+pow((latimin+j*dgr-slati[ista])*dislati,2))+1e-20;
            grdt[ista]+=gradtr[i][j]/dis;
            azi[ista]+=azig[i][j]/dis;
            weight+=1/dis;
         }
      if(ntmpp<2) { azi[ista]=-1; continue; }
      grdt[ista]=grdt[ista]/weight;
      azi[ista]=azi[ista]/weight;
   }
   return 1;
}

int calc_grdt_amp(char *file, double *slong, double *slati, double *sdat, double *grdt, double *azi, int *nsta)
{
   FILE *ff;
   char buff[300], buff2[300];
   float longmin=360, longmax=0, latimin=90, latimax=-90, longtmp, latitmp;
   double radius=6371.1391285, pi=4.0*atan(1.0), dgr=0.05, dattmp;
   double grdtx, grdty, dis, alpha, weight;
   int i, j, ii, jj, ista, nrand, ntmpp;

   if((ff=fopen(file,"r"))==NULL){
     printf("Can't open %s to read!\n",file);
     return 0;
   }
   for(*nsta=0;;){
      if(fgets(buff, 300, ff) == NULL ) break;
      if((sscanf(buff,"%lf %lf %lf", &slong[*nsta], &slati[*nsta], &sdat[*nsta]))!=3) {
         cout<<"Wrong format! Stopped: "<<buff<<endl;
         break;
      }
      if(longmin>slong[*nsta]) longmin=slong[*nsta];
      if(longmax<slong[*nsta]) longmax=slong[*nsta];
      if(latimin>slati[*nsta]) latimin=slati[*nsta];
      if(latimax<slati[*nsta]) latimax=slati[*nsta];
      *nsta=*nsta+1;
   }
   fclose(ff);
   longmin=floor(longmin/dgr)*dgr;
   latimin=floor(latimin/dgr)*dgr;
   longmax=ceil(longmax/dgr)*dgr;
   latimax=ceil(latimax/dgr)*dgr;
   int npts_long=int((longmax-longmin)/dgr+1);
   int npts_lati=int((latimax-latimin)/dgr+1);
   double dislong[npts_lati], dislati, distmp;
   double gradtr[npts_long][npts_lati], grdtmp, azig[npts_long][npts_lati], dat[npts_long][npts_lati];
   for(i=0;i<npts_long;i++) for(j=0;j<npts_lati;j++){
      gradtr[i][j] = 1e10;
      dat[i][j] = -1;
   }
   for(i=1;i<npts_lati-1;i++){
      dislati=atan(0.993277*tan((latimin+i*dgr)/180*pi))*180/pi;
      dislong[i]=radius*sin((90-dislati)/180*pi)*dgr/180*pi;
   }
   dislati=radius*dgr/180*pi;

   srand ( time(NULL) );
   nrand=rand();

   for(ista=0;ista<*nsta;ista++){
      sprintf(buff,"region_%d",nrand);
      ff=fopen(buff,"w");
      fprintf(ff, "-R%d/%d/%d/%d\0", (int)floor(longmin), (int)ceil(longmax), (int)floor(latimin), (int)ceil(latimax));
      fclose(ff);
      sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel_positive %s region_%d %f", file, nrand, dgr);
      system(buff);

      sprintf(buff,"%s.HD",file);
      if((ff=fopen(buff,"r"))==NULL){
         printf("Can't open %s to read!\n",buff);
         return 0;
      }
      for(;;){
         if( fgets(buff, 300, ff) == NULL ) break;
         sscanf(buff,"%f %f %lf", &longtmp, &latitmp, &dattmp);
         i=int((longtmp-longmin)/dgr+0.1);
         j=int((latitmp-latimin)/dgr+0.1);
         if(i<0||i>=npts_long||j<0||j>=npts_lati) continue;
         dat[i][j]=dattmp;
      }
      fclose(ff);

//   for(ista=0;ista<*nsta;ista++){
      ii=int((slong[ista]-longmin)/dgr);
      jj=int((slati[ista]-latimin)/dgr);
      for(i=ii;i<ii+2;i++)
         for(j=jj;j<jj+2;j++){
            if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
            if(gradtr[i][j] != 1e10) continue;
            grdtx=(dat[i+1][j]-dat[i-1][j])/2.0/dislong[j];
            grdty=(dat[i][j+1]-dat[i][j-1])/2.0/dislati;
            gradtr[i][j]=sqrt(grdtx*grdtx+grdty*grdty);
            azig[i][j]=pi/2-atan2(grdty,grdtx);
            if(azig[i][j]<0-1e-8)azig[i][j]+=2*pi;
            distmp=sqrt(dislong[j]*dislong[j]+dislati*dislati);
            grdtx=(dat[i+1][j-1]-dat[i-1][j+1])/2.0/distmp;
            grdty=(dat[i+1][j+1]-dat[i-1][j-1])/2.0/distmp;
            alpha=2*atan(dislati/dislong[j]);
            alpha=atan((grdtx/grdty-cos(alpha))/sin(alpha));
            grdtmp=fabs(grdty/cos(alpha));
            //cout<<gradtr[i][j]<<" "<<grdtmp<<endl;
            if(fabs((gradtr[i][j]-grdtmp)/gradtr[i][j])>0.15) gradtr[i][j]=1e20;
            else gradtr[i][j]=(gradtr[i][j]+grdtmp)/2.;
         }
   }

   for(ista=0;ista<*nsta;ista++){
      ii=int((slong[ista]-longmin)/dgr);
      jj=int((slati[ista]-latimin)/dgr);
      ntmpp=0; weight=0; grdt[ista]=0; azi[ista]=0;
      for(i=ii;i<ii+2;i++)
         for(j=jj;j<jj+2;j++){
            if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
            if(gradtr[i][j]>=1e10) continue;
            ntmpp++;
            dis=sqrt(pow((longmin+i*dgr-slong[ista])*dislong[ista],2)+pow((latimin+j*dgr-slati[ista])*dislati,2))+1e-20;
            grdt[ista]+=gradtr[i][j]/dis;
            azi[ista]+=azig[i][j]/dis;
            weight+=1/dis;
         }
      if(ntmpp<2) { azi[ista]=-1; continue; }
      grdt[ista]=grdt[ista]/weight;
      azi[ista]=azi[ista]/weight;
   }
   return 1;
}


int calc_lplc(char *file, double *slong, double *slati, double *grdt, double *azi, double *lplc, int *nsta)
{
   FILE *ff;
   char buff[300], buff2[300];
   float longmin=360, longmax=0, latimin=90, latimax=-90, longtmp, latitmp;
   double radius=6371.1391285, pi=4.0*atan(1.0), dgr=0.05, dattmp;
   double grdtx, grdty, dis, alpha;
   double sdat[2000], weight;
   int i, j, ii, jj, ista, nrand, ntmpp;

   if((ff=fopen(file,"r"))==NULL){
     printf("Can't open %s to read!\n",file);
     return 0;
   }
   for(*nsta=0;;){
      if(fgets(buff, 300, ff) == NULL ) break;
      if((sscanf(buff,"%lf %lf %lf", &slong[*nsta], &slati[*nsta], &sdat[*nsta]))!=3) {
         cout<<"Wrong format! Stopped: "<<buff<<endl;
         break;
      }
      if(longmin>slong[*nsta]) longmin=slong[*nsta];
      if(longmax<slong[*nsta]) longmax=slong[*nsta];
      if(latimin>slati[*nsta]) latimin=slati[*nsta];
      if(latimax<slati[*nsta]) latimax=slati[*nsta];
      *nsta=*nsta+1;
   }
   fclose(ff);
   longmin=floor(longmin/dgr)*dgr;
   latimin=floor(latimin/dgr)*dgr;
   longmax=ceil(longmax/dgr)*dgr;
   latimax=ceil(latimax/dgr)*dgr;
   int npts_long=int((longmax-longmin)/dgr+1);
   int npts_lati=int((latimax-latimin)/dgr+1);
   double dislong[npts_lati], dislati, distmp;
   double gradtr[npts_long][npts_lati], grdtmp, gradx[npts_long][npts_lati], grady[npts_long][npts_lati], grdx[*nsta], grdy[*nsta], azig[npts_long][npts_lati], dat[npts_long][npts_lati];

   for(i=0;i<npts_long;i++) for(j=0;j<npts_lati;j++){
      gradtr[i][j] = 1e10;
      dat[i][j] = -1;
   }
   for(i=1;i<npts_lati-1;i++){
      dislati=atan(0.993277*tan((latimin+i*dgr)/180*pi))*180/pi;
      dislong[i]=radius*sin((90-dislati)/180*pi)*dgr/180*pi;
   }
   dislati=radius*dgr/180*pi;

   nrand=rand();
   sprintf(buff,"region_%d",nrand);
   ff=fopen(buff,"w");
   fprintf(ff, "-R%d/%d/%d/%d\0", (int)floor(longmin), (int)ceil(longmax), (int)floor(latimin), (int)ceil(latimax));
   fclose(ff);
   sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel_positive %s region_%d %f", file, nrand, dgr);
   system(buff);

   sprintf(buff,"%s.HD",file);
   if((ff=fopen(buff,"r"))==NULL){
      printf("Can't open %s to read!\n",buff);
      return 0;
   }
   for(;;){
      if( fgets(buff, 300, ff) == NULL ) break;
      sscanf(buff,"%f %f %lf", &longtmp, &latitmp, &dattmp);
      i=int((longtmp-longmin)/dgr+0.1);
      j=int((latitmp-latimin)/dgr+0.1);
      if(i<0||i>=npts_long||j<0||j>=npts_lati) continue;
      dat[i][j]=dattmp;
   }
   fclose(ff);
   for(ista=0;ista<*nsta;ista++){
      ii=int((slong[ista]-longmin)/dgr);
      jj=int((slati[ista]-latimin)/dgr);
      for(i=ii;i<ii+2;i++)
         for(j=jj;j<jj+2;j++){
            if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
            if(gradtr[i][j] != 1e10) continue;
            gradx[i][j]=(dat[i+1][j]-dat[i-1][j])/2.0/dislong[j];
            grady[i][j]=(dat[i][j+1]-dat[i][j-1])/2.0/dislati;
            gradtr[i][j]=sqrt(pow(gradx[i][j],2)+pow(grady[i][j],2));
            azig[i][j]=pi/2-atan2(grady[i][j],gradx[i][j]);
            if(azig[i][j]<0-1e-8)azig[i][j]+=2*pi;
            distmp=sqrt(dislong[j]*dislong[j]+dislati*dislati);
            grdtx=(dat[i+1][j-1]-dat[i-1][j+1])/2.0/distmp;
            grdty=(dat[i+1][j+1]-dat[i-1][j-1])/2.0/distmp;
            alpha=2*atan(dislati/dislong[j]);
            alpha=atan((grdtx/grdty-cos(alpha))/sin(alpha));
            grdtmp=fabs(grdty/cos(alpha));
            if(fabs((gradtr[i][j]-grdtmp)/gradtr[i][j])>0.15) gradtr[i][j]=1e20;
            else gradtr[i][j]=(gradtr[i][j]+grdtmp)/2.;
         }
   }

   for(ista=0;ista<*nsta;ista++){
      ii=int((slong[ista]-longmin)/dgr);
      jj=int((slati[ista]-latimin)/dgr);
      ntmpp=0; weight=0; grdx[ista]=0, grdy[ista]=0, grdt[ista]=0; azi[ista]=0;
      for(i=ii;i<ii+2;i++)
         for(j=jj;j<jj+2;j++){
            if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
            if(gradtr[i][j]>=1e10) continue;
            ntmpp++;
            dis=sqrt(pow((longmin+i*dgr-slong[ista])*dislong[ista],2)+pow((latimin+j*dgr-slati[ista])*dislati,2))+1e-20;
            grdx[ista]+=gradx[i][j]/dis;
            grdy[ista]+=grady[i][j]/dis;
            grdt[ista]+=gradtr[i][j]/dis;
            azi[ista]+=azig[i][j]/dis;
            weight+=1/dis;
         }
      if(ntmpp<2){ grdt[ista]=1e10; continue; }
      grdx[ista]=grdx[ista]/weight;
      grdy[ista]=grdy[ista]/weight;
      grdt[ista]=grdt[ista]/weight;
      azi[ista]=azi[ista]/weight;
   }

   sprintf(buff,"%s_grdx",file);
   if((ff=fopen(buff,"w"))==NULL){
      printf("Can't open %s to write!\n",buff);
      return 0;
   }
   for(ista=0;ista<*nsta;ista++)
      if(grdt[ista]!=1e10) fprintf(ff,"%8.4f  %8.4f  %8g\n",slong[ista], slati[ista], grdx[ista]);
   fclose(ff);

   sprintf(buff,"%s_grdy",file);
   if((ff=fopen(buff,"w"))==NULL){
      printf("Can't open %s to write!\n",buff);
      return 0;
   }
   for(ista=0;ista<*nsta;ista++)
      if(grdt[ista]!=1e10) fprintf(ff,"%8.4f  %8.4f  %8g\n",slong[ista], slati[ista], grdy[ista]);
   fclose(ff);

   sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel %s_grdx region_%d %f", file, nrand, dgr);
   system(buff);
   sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel %s_grdy region_%d %f", file, nrand, dgr);
   system(buff);

   for(i=0;i<npts_long;i++) for(j=0;j<npts_lati;j++){
      gradtr[i][j] = 1e10;
      gradx[i][j] = -1;
   }
   
   sprintf(buff,"%s_grdx.HD",file);
   if((ff=fopen(buff,"r"))==NULL){
      printf("Can't open %s to read!\n",buff);
      return 0;
   }
   for(;;){
      if( fgets(buff, 300, ff) == NULL ) break;
      sscanf(buff,"%f %f %lf", &longtmp, &latitmp, &dattmp);
      i=int((longtmp-longmin)/dgr+0.1);
      j=int((latitmp-latimin)/dgr+0.1);
      if(i<0||i>=npts_long||j<0||j>=npts_lati) continue;
      gradx[i][j]=dattmp;
   }
   fclose(ff);

   sprintf(buff,"%s_grdy.HD",file);
   if((ff=fopen(buff,"r"))==NULL){
      printf("Can't open %s to read!\n",buff);
      return 0;
   }
   for(;;){
      if( fgets(buff, 300, ff) == NULL ) break;
      sscanf(buff,"%f %f %lf", &longtmp, &latitmp, &dattmp);
      i=int((longtmp-longmin)/dgr+0.1);
      j=int((latitmp-latimin)/dgr+0.1);
      if(i<0||i>=npts_long||j<0||j>=npts_lati) continue;
      grady[i][j]=dattmp;
   }
   fclose(ff);

   for(ista=0;ista<*nsta;ista++){
      ii=int((slong[ista]-longmin)/dgr);
      jj=int((slati[ista]-latimin)/dgr);
      for(i=ii;i<ii+2;i++)
         for(j=jj;j<jj+2;j++){
            if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
            if(gradtr[i][j] != 1e10) continue;
            grdtx=(gradx[i+1][j]-gradx[i-1][j])/2.0/dislong[j];
            grdty=(grady[i][j+1]-grady[i][j-1])/2.0/dislati;
            grdtmp=grdtx+grdty;
            alpha=atan(dislong[j]/dislati);
            gradtr[i][j]=-gradx[i-1][j+1]*sin(alpha)+grady[i-1][j+1]*cos(alpha)+grady[i][j+1]+gradx[i+1][j+1]*sin(alpha)+grady[i+1][j+1]*cos(alpha)-gradx[i-1][j]+gradx[i+1][j]-gradx[i-1][j-1]*sin(alpha)-grady[i-1][j-1]*cos(alpha)-grady[i][j-1]+gradx[i+1][j-1]*sin(alpha)-grady[i+1][j-1]*cos(alpha);
            gradtr[i][j]=gradtr[i][j]/8*(dislati+dislong[j])/dislati/dislong[j];
            if(fabs((grdtmp-gradtr[i][j])/gradtr[i][j])>0.3) gradtr[i][j]=1e20;
            //else gradtr[i][j]=grdtmp;
         }
   }

   for(ista=0;ista<*nsta;ista++){
      ii=int((slong[ista]-longmin)/dgr);
      jj=int((slati[ista]-latimin)/dgr);
      ntmpp=0; weight=0; lplc[ista]=0;
      for(i=ii;i<ii+2;i++)
         for(j=jj;j<jj+2;j++){
            if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
            if(gradtr[i][j]>=1e10) continue;
            ntmpp++;
            dis=sqrt(pow((longmin+i*dgr-slong[ista])*dislong[ista],2)+pow((latimin+j*dgr-slati[ista])*dislati,2))+1e-20;
            lplc[ista]+=gradtr[i][j]/dis;
            weight+=1/dis;
         }
      if(ntmpp<2) { azi[ista]=-1; continue; }
      lplc[ista]=lplc[ista]/weight;
   }

   sprintf(buff,"rm -f %s_grdx %s_grdx.HD %s_grdy %s_grdy.HD region_%d", file, file, file, file, nrand);
   system(buff);

   return 1;
}

int main (int argc, char *argv[])
{
   if(argc!=3){
      cout<<"usage:csta_map_gradient [input_trvt_map] [input_amp_map]"<<endl;
      return 0;
   }
   FILE *ff;
   char fname[100], *tmpc;
   double grdt[2000], azit[2000], lplc[2000], amp[2000], grda[2000], azia[2000];
   double tlong[2000], tlati[2000], along[2000], alati[2000];
   int i, nsta, nstt;
   calc_grdt_amp(argv[2], along, alati, amp, grda, azia, &nsta);
   calc_lplc(argv[1], tlong, tlati, grdt, azit, lplc, &nstt);
   if(nsta!=nstt){
      cout<<"Different station numbers from phase and amplitude map!"<<endl;
      return 0;
   }
   double dec[nsta];
   for(i=0;i<nsta;i++){
      if(azia[i]==-1||azit[i]==-1) continue;
      dec[i]=2*grda[i]*grdt[i]*cos(azia[i]-azit[i])/amp[i];
   }

   tmpc = strdup(argv[2]);
   tmpc = strtok(tmpc, "_");

   sprintf(fname,"%s_apr_amp_decay",tmpc);
   ff=fopen(fname,"w");
   for(i=0;i<nsta;i++){
      if(azia[i]==-1||azit[i]==-1) continue;
      fprintf(ff,"%8.4f  %8.4f  %8g\n", along[i], alati[i], dec[i]);
   }

   sprintf(fname,"%s_fcs_dfcs",tmpc);
   ff=fopen(fname,"w");
   for(i=0;i<nsta;i++){
      if(azia[i]==-1||azit[i]==-1) continue;
      fprintf(ff,"%8.4f  %8.4f  %8g\n", along[i], alati[i], lplc[i]);
   }
/*
   sprintf(fname,"%s_tgrdt",tmpc);
   ff=fopen(fname,"w");
   for(i=0;i<nsta;i++){
      if(azia[i]==-1||azit[i]==-1) continue;
      fprintf(ff,"%8.4f  %8.4f  %8g\n", along[i], alati[i], grdt[i]);
   }
*/
   sprintf(fname,"%s_atn_aplf_term",tmpc);
   ff=fopen(fname,"w");
   for(i=0;i<nsta;i++){
      if(azia[i]==-1||azit[i]==-1) continue;
      fprintf(ff,"%8.4f  %8.4f  %8g\n", along[i], alati[i], dec[i]+lplc[i]);
   }

   return 1;
}
