#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

main(int na, char *arg[])
{
  if (na!=6)
    {
      cout<<"usage:get_dist lati_1 long_1 lati_2 long_2 [d/a1/a2/b]"<<endl;
      return 0;
    }
  int i,flag=0;
  double lati1,long1,lati2,long2,latio1,longo1,dlt_lati,dlt_long;
  double Ra,Rb,R,Rx,f;
  double U1,U2;
  double s,Ds,alpha1,alpha2,theta;
  double ctr_angl,dist;
  double pi;
  double cv,cv1,cv2,cv3,cv4,cv5,cvC,numda1;
  double numda;
  double mius, cvA, cvB, deltacv;

  pi=4.0*atan(1.0);
  Ra = 6378.137; Rb = 6356.7523142;
  f = 1/298.257223563;
  lati1=atof(arg[1]);
  lati2=atof(arg[3]);
  long1=atof(arg[2]);
  long2=atof(arg[4]);
begin:
  long1=long1-(int)floor(long1)/360*360;
  if(long1<0) long1+=360;
  long2=long2-(int)floor(long2)/360*360;
  if(long2<0) long2+=360;
  if(lati1==-lati2 && fabs(long1-long2)==180){
/*
  latio1=lati1*pi/180.;
  Rx=sqrt((pow(Ra*Ra*cos(latio1),2)+pow(Rb*Rb*sin(latio1),2))/(pow(Ra*cos(latio1),2)+pow(Rb*sin(latio1),2)));
  s=pi*(Rx+Ra)/2.;
*/
  s=20003.917357;
  alpha1=999;
  alpha2=999;
  goto output;
  }
  dlt_lati=fabs(lati2-lati1);
  dlt_long=long2-long1;

if (dlt_long > 180.000)  dlt_long = 360.000000 - dlt_long;
if (dlt_long < -180.000) dlt_long = 360.000 - fabs(dlt_long);
dlt_long = fabs(dlt_long);

U1 = atan((1-f)*tan(lati1/180*pi));
U2 = atan((1-f)*tan(lati2/180*pi));
dlt_long = dlt_long*pi/180;

if(flag==0){
alpha1=atan2(cos(U2)*sin(dlt_long), cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(dlt_long))*180/pi;
alpha2=atan2(cos(U1)*sin(dlt_long), -sin(U1)*cos(U2) + cos(U1)*sin(U2)*cos(dlt_long))*180/pi;
if( fabs(long2-long1)>180 ) { alpha1 = 360-alpha1; alpha2 = 360-alpha2;}
if( long2 < long1 ) { alpha1 = 360-alpha1; alpha2 = 360-alpha2; }
}

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
//printf("cv1: %lf  cv2: %lf  numda: %lf  cv: %lf\n",cv1,cv2,numda,cv);
cv3 = cos(U1)*cos(U2)*sin(numda)/sin(cv);   
//printf("cv3: %f\n",cv3);
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
//printf("%lf\n",Rx);
goto begin;
}
} while (fabs(numda - numda1) > 1e-10);
//printf("cv4: %f\n",cv4);
mius = cv4*(Ra*Ra - Rb*Rb)/(Rb*Rb);
cvA = 1+mius/16384*(4096 + mius*(-768 + mius*(320 - 175*mius)));
//cout<<"cvA "<<cvA<<endl;
cvB = mius/1024*(256+ mius*(-128 + mius*(74 - 47*mius)));
//printf("cvB: %f  cv1: %f  cv2: %f  cv5: %f\n",cvB,cv1,cv2,cv5);
deltacv = cvB*cv1*(cv5 +cvB/4*(cv2*(-1 + 2*cv5*cv5)-cvB/6*cv5*(-3+4*cv1*cv1)*(-3+4*cv5*cv5) ));
//cout<<"delatacv "<<deltacv<<"cvb "<<cvb<<"cv "<<cv<<endl;
//printf("cvA:%f deltacv:%f mius:%f\n",cvA,deltacv,mius);
s = Rb * cvA *(cv - deltacv);

if(flag==1){
theta=alpha1;
if(theta>180) theta=360-theta;
theta=fabs(90-theta);
Ds=Rx*(90-theta)/90.+(Ra+Rb)*pi*theta/180.;
s=Ds-s;
}

output:
if(strcmp(arg[5],"b")==0) printf("%lf\na1: %lf  a2: %lf\n",s,alpha1,alpha2);
else if(strcmp(arg[5],"d")==0) printf ("%lf\n",s);
else if(strcmp(arg[5],"a1")==0) printf("%lf\n",alpha1);
else if(strcmp(arg[5],"a2")==0) printf("%lf\n",alpha2);
else cout<<"arg[5]: enter \'b\' for both, \'d\' for dist, \'a1\' for azi1, \'a2\' for azi2"<<endl;
return 1;
}
