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
  int i, j, k, l, l2, step, slong, slati, nlong, nsta, nesta;
  char stalst[5000][10], filename[100], buff[100];
  float jf, longlst[5000], latilst[5000], longl, longh, latil, latih;
  double amptmp, amp1[5000], amp_tmp[2000], amp_azimuth[5000][3601];
  double dist[5000], distmin, alpha;
  double sp_n[360][181]={0};

  if(argc != 3)
    {
      printf("Usage: Search_ASN_Source [region_infile] [sta_lst]\n");
      exit(-1);
    }

  if((fin = fopen(argv[1],"r"))==NULL) {
    printf("Cannot open region file %s\n",argv[1]);
    exit(1);
  }
  if((fscanf(fin,"%s", buff))!=1) {
    printf("Wrong region file format!\n");
    exit(1);
  }
  sscanf(buff,"-R%f/%f/%f/%f",&longl,&longh,&latil,&latih);
  fclose(fin);

  if((fin = fopen(argv[2], "r"))==NULL) {
    printf("Cannot open station file %s\n",argv[2]);
    exit(1);
  }
  for(i=0;;i++)
    if((fscanf(fin,"%s %f %f", stalst[i], &longlst[i], &latilst[i]))!=3)
       break;
  nsta = i;
  fclose(fin);
  cout<<"nsta: "<<nsta<<endl;

  for(i=0;i<nsta;i++) for(j=0;j<3601;j++) amp_azimuth[i][j]=-1.0;
  for(i=0;i<nsta;i++){
     sprintf(filename,"%s_amp_azimuth",stalst[i]);
     if((fin=fopen(filename,"r"))==NULL) {
        cout<<"Can not find amp file "<<filename<<". Skipped"<<endl;
        continue;
     }
     amp_azimuth[i][3600]=1;
     for(;;) {
        if((fscanf(fin,"%f %lf",&jf,&amptmp))!=2) break;
        amp_azimuth[i][(int)(jf*10)] = amptmp;
     }
     fclose(fin);
  }

  nlong=(int)floor(longh)-(int)ceil(longl)+1;
  for(slong=(int)ceil(longl);slong<=(int)floor(longh);slong++) {
     i=slong-(int)ceil(longl);
//     for(slati=(int)ceil(latil)+5;slati<=(int)ceil(latil)+5;slati++) {
     for(slati=(int)ceil(latil);slati<=(int)floor(latih);slati++) {
        j=slati-(int)ceil(latil);
        sp_n[i][j]=0;
        for(k=0;k<nsta;k++) {
           if(amp_azimuth[k][3600]==-1) {amp1[k]=-1.0; continue; }
           dist_azimuth(slati, slong, latilst[k], longlst[k], &dist[k], &alpha);
           l=(int)(alpha*10+0.5);
//           if(l<alpha*10) step=1;
//           else step=-1;
//           l2=l+step;
           if(l==3600) l=0;
//           if(l2==3600) l2=0;
           amp1[k]=amp_azimuth[k][l];
//           if(amp1[k]==-1) amp1[k]=amp_azimuth[k][l2];
/*           l=(int)floor(alpha+0.5);
           if(l<alpha) step=1;
           else step=-1;
           l2=l+step;
           if(l==360) l=0;
           else if(l2==360) l2=0;
          // else if(l2==-1) l2=359;
           amp1=amp_azimuth[k][l];
           if(amp1==-1) continue;
           amp2=amp_azimuth[k][l2];
           if(amp2==-1) { 
              l2+=step;
              if(l2==360) l2=0;
              else if(l2==-1) l2=359;
              amp2=amp_azimuth[k][l2];
           }
           if(amp2!=-1) amp1=(amp1*(l2-alpha)+amp2*(alpha-l))/(l2-l);
           //sp_n[i][j][0]*=(amp1*1e6); */
//           sp_n[i][j][0]+=amp1;
//           sp_n[i][j][1]+=1; 
//cout<<"amp: "<<sp_n[i][j][0]<<"  num: "<<sp_n[i][j][1]<<endl;
        }
        nesta=0;
        for(k=0;k<nsta;k++) if(amp1[k]!=-1) nesta+=1;
        if(nesta<nsta*0.5) continue;
        for(k=0;k<nesta/1.0;k++) {
           distmin=50000.0;
           for(l=0;l<nsta;l++) {
              if(amp1[l]>-1 && dist[l]<distmin) {
                 distmin=dist[l];
                 l2=l;
              }
           }
          // if(amp1[l2]<0) cout<<amp1[l2]<<endl;
           amp_tmp[k]=amp1[l2]*amp1[l2];
           amp1[l2]=-1.0;
//        cout.precision(10);
//        cout<<"slong: "<<slong<<"  slati: "<<slati<<"  sta: "<<stalst[l2]<<"  amp: "<<sqrt(amp_tmp[k])<<endl;
        }
        nesta=k;
        for(k=0;k<nesta/50.0;k++) {
           amptmp=0.0;
           for(l=0;l<nesta;l++) {
              if( amp_tmp[l]>amptmp ) {
                 amptmp=amp_tmp[l];
                 l2=l;
              }
           }
           amp_tmp[l2]=0;
        }
        l2=nesta-k;
        for(k=0;k<nesta;k++) sp_n[i][j]+=amp_tmp[k];
        sp_n[i][j]/=(float)l2;
        //sp_n[i][j][0]=pow(sp_n[i][j][0],1./sp_n[i][j][1]);
     }
     cout.precision(1);
     cout<<fixed<<(i+1)*100./nlong<<" percent completed."<<endl;
  }
 
  sprintf(filename,"source_p_map.HD");
  if((fin = fopen(filename,"w")) ==NULL) {
     printf("Cannot open file %s to write\n",filename);
     return 0;
  }
  for(slong=(int)ceil(longl);slong<=(int)floor(longh);slong++) {
     i=slong-(int)ceil(longl);
     for(slati=(int)ceil(latil);slati<=(int)floor(latih);slati++) {
        j=slati-(int)ceil(latil);
        if(sp_n[i][j]!=0)
          fprintf(fin,"%3d  %2d  %.4g\n", slong, slati, sp_n[i][j]);
     }
  }
  fclose(fin);

}
