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


int main (int argc, char *argv[])
{
  if(argc!=5){
      cout<<"usage:csta_map_gradient [input_map] [min_dis] [csta_long] [csta_lati]"<<endl;
      return 0;
  }
  FILE *ff;
  char buff[300];
  double longmin=360, longmax=0, latimin=90, latimax=-90, longtmp, latitmp;
  double mdis=atof(argv[2]), clong=atof(argv[3]), clati=atof(argv[4]);
  double radius=6371.1391285, pi=4.0*atan(1.0), dgr=0.05, dattmp;
  double grdtx, grdty, grdavg, azi, alpha;
  double stalong[2000],stalati[2000],stadat[2000],stadis[2000];
  int i, j, ii, jj, ista, nsta, ntmp, ntmpp;
  
  if(clong<0) clong+=360;

  if((ff=fopen(argv[1],"r"))==NULL){
     printf("Can't open %s to read!\n",argv[1]);
     return 0;
  }
  for(nsta=0;;){
     if(fgets(buff, 300, ff) == NULL ) break;
     if((sscanf(buff,"%lf %lf %lf", &stalong[nsta], &stalati[nsta], &stadat[nsta]))!=3) {
        cout<<"Wrong format! Stopped: "<<buff<<endl;
        break;
     }
     calc_dist( clati, clong, stalati[nsta], stalong[nsta], &stadis[nsta] );
     if(stadis[nsta]<mdis-150) continue;
     if(stalong[nsta]<0) stalong[nsta]+=360;
     if(longmin>stalong[nsta]) longmin=stalong[nsta];
     if(longmax<stalong[nsta]) longmax=stalong[nsta];
     if(latimin>stalati[nsta]) latimin=stalati[nsta];
     if(latimax<stalati[nsta]) latimax=stalati[nsta];
     nsta++;
  }
  fclose(ff);
  longmin=floor(longmin/dgr)*dgr;
  latimin=floor(latimin/dgr)*dgr;
  longmax=ceil(longmax/dgr)*dgr;
  latimax=ceil(latimax/dgr)*dgr;
  int npts_long=int((longmax-longmin)/dgr+1);
  int npts_lati=int((latimax-latimin)/dgr+1);
  double dislong[npts_lati], dislati, distmp;
  double gradtr[npts_long][npts_lati][2], dat[npts_long][npts_lati];
//  memset(dat,-1,npts_long*npts_lati*sizeof(double));
//  memset(gradr,1e10,npts_long*npts_lati*sizeof(double));
  for(i=0;i<npts_long;i++) for(j=0;j<npts_lati;j++){
     gradtr[i][j][0] = 1e10;
     dat[i][j] = -1;
  }
  for(i=1;i<npts_lati-1;i++){
     dislati=atan(0.993277*tan((latimin+i*dgr)/180*pi))*180/pi;
     dislong[i]=radius*sin((90-dislati)/180*pi)*dgr/180*pi;
  }
  dislati=radius*dgr/180*pi;

//  if(access(buff, F_OK) != 0) {
     ff=fopen("region_tmp","w");
     fprintf(ff, "-R%d/%d/%d/%d\0", (int)floor(longmin), (int)ceil(longmax), (int)floor(latimin), (int)ceil(latimax));
     fclose(ff);
     sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel_positive %s region_tmp %f", argv[1], dgr);
     system(buff);
//  }

  sprintf(buff,"%s.HD",argv[1]);
  if((ff=fopen(buff,"r"))==NULL){
     printf("Can't open %s to read!\n",buff);
     return 0;
  }
  for(;;){
     if( fgets(buff, 300, ff) == NULL ) break;
     sscanf(buff,"%lf %lf %lf", &longtmp, &latitmp, &dattmp);
     i=int((longtmp-longmin)/dgr+0.1);
     j=int((latitmp-latimin)/dgr+0.1);
     if(i<0||i>=npts_long||j<0||j>=npts_lati) continue;
     dat[i][j]=dattmp;
  }
  fclose(ff);

//  grdavg=0; ntmp=0;
  for(ista=0;ista<nsta;ista++){
     if(stadis[ista]<mdis) continue;
     ii=int((stalong[ista]-longmin)/dgr);
     jj=int((stalati[ista]-latimin)/dgr);
     for(i=ii;i<ii+2;i++)
        for(j=jj;j<jj+2;j++){
           if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
           if(gradtr[i][j][0] != 1e10) continue;
           calc_azimuth(clati, clong, latimin+j*dgr, longmin+i*dgr, &azi);
           grdtx=(dat[i+1][j]-dat[i-1][j])/2.0/dislong[j];
           grdty=(dat[i][j+1]-dat[i][j-1])/2.0/dislati;
           gradtr[i][j][0]=grdtx*sin(azi*pi/180.)+grdty*cos(azi*pi/180.);
           gradtr[i][j][1]=sqrt(grdtx*grdtx+grdty*grdty);
           distmp=sqrt(dislong[j]*dislong[j]+dislati*dislati);
           grdtx=(dat[i+1][j-1]-dat[i-1][j+1])/2.0/distmp;
           grdty=(dat[i+1][j+1]-dat[i-1][j-1])/2.0/distmp;
           alpha=2*atan(dislati/dislong[j]);
           alpha=atan((grdtx/grdty-cos(alpha))/sin(alpha));
           gradtr[i][j][1]=(gradtr[i][j][1]-fabs(grdty/cos(alpha)))/gradtr[i][j][1];
           if(fabs(gradtr[i][j][1])>0.15) gradtr[i][j][0]=1e20;
           
           //grdavg+=fabs(gradtr[i][j][0]);
           //ntmp++;
        }
  }
//  grdavg/=ntmp;
  
  sprintf(buff,"%s_grdt",argv[1]);
  if((ff=fopen(buff,"w"))==NULL){
     printf("Can't open %s to write!\n",buff);
     return 0;
  }
  for(ista=0;ista<nsta;ista++){
     if(stadis[ista]<mdis) continue;
     ii=int((stalong[ista]-longmin)/dgr);
     jj=int((stalati[ista]-latimin)/dgr);
     ntmpp=0;
     for(i=ii;i<ii+2;i++)
        for(j=jj;j<jj+2;j++){
           if(gradtr[i][j][0]>=1e10) continue;
           ntmpp++;
        }
     if(ntmpp<3) {cout<<ntmpp<<endl; continue; }
     for(i=ii;i<ii+2;i++)
        for(j=jj;j<jj+2;j++){
           if(gradtr[i][j][0]>=1e10) continue;
           fprintf(ff,"%6.2f  %6.2f  %8g\n",longmin+i*dgr, latimin+j*dgr, gradtr[i][j][0]);
     }
  }
  fclose(ff);

/*
  sprintf(buff,"%s_grdt",argv[1]);
  if((ff=fopen(buff,"w"))==NULL){
     printf("Can't open %s to write!\n",buff);
     return 0;
  }
  for(i=1;i<npts_long-1;i++)
     for(j=1;j<npts_lati-1;j++) {
        if(gradtr[i][j][0]==1e10) continue;
        if(fabs((gradtr[i][j][0]-gradtr[i][j][1])/gradtr[i][j][0])<0.2 || fabs(gradtr[i][j][0]-gradtr[i][j][1])<grdavg/20.) fprintf(ff,"%6.2f  %6.2f  %8g\n",longmin+i*dgr, latimin+j*dgr, (gradtr[i][j][0]+gradtr[i][j][1])/2.);
     }
  fclose(ff);
*/
  return 1;
}
