#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
using namespace std;

int dist_azimuth(double lati1, double long1, double lati2, double long2, double *dist, double *alpha2)
{
  double dlt_lati,dlt_long;
  double Ra,Rb,R,f;
  double U1,U2;
  double ctr_angl;
  double pi;
  pi=4.0*atan(1.0);
  Ra = 6378.137; Rb = 6356.7523142;
  f = 1/298.257223563;
//  lati1=atof(arg[1]);
//  lati2=atof(arg[3]);
//  long1=atof(arg[2]);
//  long2=atof(arg[4]);
//  avg_lati=(lati1+lati2)/2;
  long1=long1-(int)floor(long1)/360*360;
  if(long1<0) long1+=360;
  long2=long2-(int)floor(long2)/360*360;
  if(long2<0) long2+=360;
  dlt_lati=fabs(lati2-lati1);
  dlt_long=long2-long1;
if (dlt_long > 180.000)  dlt_long = 360.000000 - dlt_long;
if (dlt_long < -180.000) dlt_long = 360.000 - fabs(dlt_long);
dlt_long = fabs(dlt_long);

U1 = atan((1-f)*tan(lati1/180*pi));
U2 = atan((1-f)*tan(lati2/180*pi));
dlt_long = dlt_long*pi/180;

double cv,cv1,cv2,cv3,cv4,cv5,cvC,numda1;
double numda = dlt_long;
numda1 = numda;
do {
numda = numda1;
cv1 =  sqrt( (cos(U2)*sin(numda))*(cos(U2)*sin(numda))+ (cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda))*(cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda)) );
cv2 = sin(U1)*sin(U2)+ cos(U1)*cos(U2)*cos(numda);
cv = atan2(cv1,cv2);
if(cv==0)cv+=1e-20;
//printf("cv1: %lf  cv2: %lf  numda: %lf  cv: %lf\n",cv1,cv2,numda,cv);
cv3 = cos(U1)*cos(U2)*sin(numda)/sin(cv);
//printf("cv3: %f\n",cv3);
cv4 = 1 - cv3*cv3;
if(cv4==0)cv4+=1e-20;
cv5 = cos(cv) - 2*sin(U1)*sin(U2)/cv4;
cvC = f/16*cv4*(4 + f*(4 - 3*cv4));
numda1 = dlt_long + (1-cvC)*f*cv3*(cv + cvC*cv1*(cv5 + cvC*cv2*(-1 +2*cv5*cv5)));
} while (fabs(numda - numda1) > 1e-10);
double mius, cvA, cvB, deltacv;
//double alpha1,alpha2;
//printf("cv4: %f\n",cv4);
mius = cv4*(Ra*Ra - Rb*Rb)/(Rb*Rb);
cvA = 1+mius/16384*(4096 + mius*(-768 + mius*(320 - 175*mius)));
//cout<<"cvA "<<cvA<<endl;
cvB = mius/1024*(256+ mius*(-128 + mius*(74 - 47*mius)));
//printf("cvB: %f  cv1: %f  cv2: %f  cv5: %f\n",cvB,cv1,cv2,cv5);
deltacv = cvB*cv1*(cv5 +cvB/4*(cv2*(-1 + 2*cv5*cv5)-cvB/6*cv5*(-3+4*cv1*cv1)*(-3+4*cv5*cv5) ));
//cout<<"delatacv "<<deltacv<<"cvb "<<cvb<<"cv "<<cv<<endl;
//printf("cvA:%f deltacv:%f mius:%f\n",cvA,deltacv,mius);
*dist = Rb * cvA *(cv - deltacv);
//*alpha1=atan2(cos(U2)*sin(dlt_long), cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(dlt_long))*180/pi;
*alpha2=atan2(cos(U1)*sin(dlt_long), -sin(U1)*cos(U2) + cos(U1)*sin(U2)*cos(dlt_long))*180/pi;
if( fabs(long2-long1)>180 ) *alpha2 = 360-*alpha2;
if( long2 < long1 ) *alpha2 = 360-*alpha2;
//printf ("%lf\n",s);
//printf("a1: %lf  a2: %lf\n",alpha1,alpha2);
return 1;
}


int main (int argc, char *argv[])
{
  FILE *fin;
  char stalst[5000][10], filename[100];
  int i, nsta;
  float jf, longlst[5000], latilst[5000], amplst[5000];
  double slong, slati, amptmp, ampl, amph, dist, alpha, azil, azih;

  if(argc != 3)
    {
      printf("Usage: Calc_source_amp [source_long] [source_lati]\n");
      exit(-1);
    }

  slong=atof(argv[1]);
  slati=atof(argv[2]);
/*
  azil=atof(argv[3]);
  azih=atof(argv[4]);
  azil=azil-(int)floor(azil)/360*360;
  if(azil<0) azil+=360;
  azih=azih-(int)floor(azih)/360*360;
  if(azih<0)azih+=360;
  if(azih<azil) { 
     aflag=1; 
     alpha=(azih+360.0+azil)/2.0;
     if(alpha>=360) alpha-=360;
  }
  else alpha=(azih+azil)/2.0;
*/
  if((fin = fopen("station.lst","r"))==NULL) {
    printf("Cannot open station.lst. Change to the Amp_Azi_per dir before running\n");
    exit(1);
  }
  for(i=0;;i++) {
    if((fscanf(fin,"%s %f %f", stalst[i], &longlst[i], &latilst[i]))!=3)
       break;
    if(longlst[i]<0) longlst[i]+=360;
  }
  nsta = i;
  fclose(fin);
  cout<<"nsta: "<<nsta<<endl;

//  sprintf(filename,"Source_Amp_Map%f_%f",slong,slati);
//  system("rm -rf %s",filename);
//  system("mkdir %s",filename);
  
  for(i=0;i<nsta;i++){
     amplst[i]=-1;
     dist_azimuth(slati, slong, latilst[i], longlst[i], &dist, &alpha);
     azil=(int)(alpha*10)/10.0;
     azih=azil+0.1;
     sprintf(filename,"%s_amp_azimuth",stalst[i]);
     if((fin=fopen(filename,"r"))==NULL) {
        cout<<"Can not find amp file "<<filename<<". Skipped"<<endl;
        continue;
     }
     ampl=0; amph=0;
     for(;;) {
        if((fscanf(fin,"%f %lf",&jf,&amptmp))!=2) break;
        if(jf<azil) continue;
        if(jf==azil) { ampl = amptmp; continue; }
        if(jf==azih) { amph = amptmp; break; }
        if(jf>azih) break;
     }
     fclose(fin);
     if(ampl==0 && amph==0) continue;
     else if(ampl==0) amplst[i]=amph;
     else if(amph==0) amplst[i]=ampl;
     else amplst[i]=(ampl*(azih-alpha)+amph*(alpha-azil))*10;
  }

  sprintf(filename,"Source_Amp_Map_%.1f_%.1f",slong,slati);
  if((fin = fopen(filename,"w")) ==NULL) {
     printf("Cannot open file %s to write\n",filename);
     return 0;
  }
  for(i=0;i<nsta;i++){
     if(amplst[i]==-1) continue;
     fprintf(fin,"%7.3f  %7.3f  %.4g\n", longlst[i], latilst[i], amplst[i]);
  }
  fclose(fin);

}
