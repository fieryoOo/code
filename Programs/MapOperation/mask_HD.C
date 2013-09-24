#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
using namespace std;

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

int main(int argc, char *argv[])
{
   if(argc!=3){
      cout<<"Usage: mask_HD [in_HD_file] [in_txt_file]"<<endl;
      return 0;
   }
   FILE *ff;
   char buff[300];
   float ib,ie, jb,je, longmin=360, longmax=0, latimin=90, latimax=-90, tmp;
   float dgr=100, stalon[200000], stalat[200000], grdlon[200000], grdlat[200000];
   double radius=6371.1391285, pi=4.0*atan(1.0), grddat[200000];
   double dislati, dislong;
   int i, j, ista, nsta, igrid, ngrid;

   if((ff=fopen(argv[2],"r")) == NULL){
      cout<<"Can't open file "<<argv[2]<<" for read!"<<endl;
      return 0;
   }
   for(nsta=0;;nsta++){
      if(fgets(buff, 300, ff) == NULL) break;
      sscanf(buff,"%f %f", &stalon[nsta], &stalat[nsta]);
      if(stalon[nsta]<0) stalon[nsta]+=360;
   }
   fclose(ff);
   if((ff=fopen(argv[1],"r")) == NULL){
      cout<<"Can't open file "<<argv[1]<<" for read!"<<endl;
      return 0;
   }
   for(ngrid=0;;ngrid++){
      if(fgets(buff, 300, ff) == NULL) break;
      sscanf(buff,"%f %f %lf", &grdlon[ngrid], &grdlat[ngrid], &grddat[ngrid]);
      if(longmin>grdlon[ngrid]) longmin=grdlon[ngrid];
      if(longmax<grdlon[ngrid]) longmax=grdlon[ngrid];
      if(latimin>grdlat[ngrid]) latimin=grdlat[ngrid];
      if(latimax<grdlat[ngrid]) latimax=grdlat[ngrid];
   }
   fclose(ff);
   if(ngrid<100){
      cout<<"Only "<<ngrid<<" points in HD file!"<<endl;
      return 0;
   }

   for(i=0;i<10;i++){
      tmp=fabs(grdlon[i]-grdlon[i+1]);
      if(tmp==0) tmp=fabs(grdlat[i]-grdlat[i+1]);
      if(tmp==0) continue;
      if(dgr>tmp) dgr=tmp;
   }
   if(dgr==100){
      cout<<"Wrong HD file!"<<endl;
      return 0;
   }

   longmin=floor(longmin/dgr)*dgr;
   latimin=floor(latimin/dgr)*dgr;
   longmax=ceil(longmax/dgr)*dgr;
   latimax=ceil(latimax/dgr)*dgr;
   int npts_long=int((longmax-longmin)/dgr+1);
   int npts_lati=int((latimax-latimin)/dgr+1);
   int flag[npts_long][npts_lati];
   double dat[npts_long][npts_lati];

   for(i=0;i<npts_long;i++) for(j=0;j<npts_lati;j++){
     flag[i][j] = -1;
     dat[i][j] = -1;
   }

   for(igrid=0;igrid<ngrid;igrid++){
      i=int((grdlon[igrid]-longmin)/dgr+0.1);
      j=int((grdlat[igrid]-latimin)/dgr+0.1);
      dat[i][j]=grddat[igrid];
   }

   float edis=50.;
   dislati=radius*dgr/180*pi;
   for(ista=0;ista<nsta;ista++){
      dislong=atan(0.993277*tan(stalat[ista]/180*pi))*180/pi;
      dislong=radius*sin((90-dislong)/180*pi)*dgr/180*pi;
      ib=(stalon[ista]-longmin)/dgr-edis/dislong;
      ie=ib+2*edis/dislong;
      if(ib<0)ib=0; if(ie>npts_long-1)ie=npts_long-1;
      jb=(stalat[ista]-latimin)/dgr-edis/dislati;
      je=jb+2*edis/dislati;
      if(jb<0)jb=0; if(je>npts_lati-1)je=npts_lati-1;
      for(i=(int)ceil(ib);i<=ie;i++)
         for(j=(int)ceil(jb);j<=je;j++)
            if(pow(((stalon[ista]-longmin)/dgr-i)*dislong,2)+pow(((stalat[ista]-latimin)/dgr-j)*dislati,2) < edis*edis) flag[i][j]=1;
   }

   sprintf(buff, "%s_mask\0", argv[1]);
   if((ff=fopen(buff,"w")) == NULL){
      cout<<"Can't open file "<<buff<<" for write!"<<endl;
      return 0;
   }
   for(i=0;i<npts_long;i++) 
      for(j=0;j<npts_lati;j++)
         if(flag[i][j]==1 && dat[i][j]!=-1)
            fprintf(ff,"%5.1f  %5.1f  %8g\n", longmin+i*dgr, latimin+j*dgr, dat[i][j]);
   fclose(ff);
   return 1;
}
