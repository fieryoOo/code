#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
using namespace std;

int dist_azimuth(double lati1, double long1, double lati2, double long2, double *dist, double *alpha1)
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
if(cv==0)cv+=0.00000000001;
//printf("cv1: %lf  cv2: %lf  numda: %lf  cv: %lf\n",cv1,cv2,numda,cv);
cv3 = cos(U1)*cos(U2)*sin(numda)/sin(cv);
//printf("cv3: %f\n",cv3);
cv4 = 1 - cv3*cv3;
if(cv4==0)cv4+=0.00000000001;
cv5 = cos(cv) - 2*sin(U1)*sin(U2)/cv4;
cvC = f/16*cv4*(4 + f*(4 - 3*cv4));
numda1 = dlt_long + (1-cvC)*f*cv3*(cv + cvC*cv1*(cv5 + cvC*cv2*(-1 +2*cv5*cv5)));
} while (fabs(numda - numda1) > 0.0000000001);
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
*alpha1=atan2(cos(U2)*sin(dlt_long), cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(dlt_long))*180/pi;
//alpha2=atan2(cos(U1)*sin(dlt_long), -sin(U1)*cos(U2) + cos(U1)*sin(U2)*cos(dlt_long))*180/pi;
if( long2 < long1 ) *alpha1 = 360-*alpha1;
//printf ("%lf\n",s);
//printf("a1: %lf  a2: %lf\n",alpha1,alpha2);
return 1;
}


int main (int argc, char *argv[])
{
  
  FILE *fin, *fout;
  char filename[300], outname[300], buff[300];
  char osta[10], stalst[5000][10], *tmp, *staname;
  char sta_pair[1000000][50], sta1[1000000][10], sta2[1000000][10];
  int i, j, k, iper, jper, nfile, nsta, deg_l, deg_h, daynum[1000000];
  float olong, olati, longlst[5000], latilst[5000];
  float per, perl, perh, pertmp[30], amptmp[30], fstep=-0.08;
  float period, amppos, snrpos, ampneg, snrneg;
  double dist, alpha, d_alpha, e_n=2.7182818284590451, powe, amp_n[360][4];
  if(argc != 4) 
    {
      printf("Usage: calc_amp_azimuth [sta.lst] [file_daynum.lst] [period]\n"); 
      exit(-1);
    }
  
  per=atof(argv[3]);
  perl=1./exp(log(1./per)-fstep);
  perh=1./exp(log(1./per)+fstep);

  if((fin = fopen(argv[2], "r"))==NULL) {
    printf("Cannot open file list %s\n",argv[2]);
    exit(1);
  }
  for(i=0;;i++) {
    if((fscanf(fin,"%s %d",sta_pair[i], &daynum[i]))!=2) break;
    tmp = strdup(sta_pair[i]);
    strtok(tmp, "_");
    staname = strtok(NULL, "_");
    strcpy(sta1[i],staname);
    staname = strtok(NULL, ".");
    strcpy(sta2[i],staname);
  }
  nfile = i;
  fclose(fin);

  if((fin = fopen(argv[1], "r"))==NULL) {
    printf("Cannot open station file %s\n",argv[1]);
    exit(1);
  }
  for(i=0;;i++) {
    if((fscanf(fin,"%s %f %f", stalst[i], &longlst[i], &latilst[i]))!=3) break;
  }
  nsta = i;
  fclose(fin);
  for(i=0;i<nsta;i++) {
//  for(i=1298;i<1299;i++) {
     memset(amp_n,0,360*4*sizeof(double));
     for(j=0;j<nfile;j++) {
        if(daynum[j]<60) continue;
        if(strcmp(stalst[i],sta1[j])==0) {
           strcpy(osta,sta2[j]);
           for(k=0;k<nsta;k++)
              if(strcmp(osta,stalst[k])==0) {
                 olong=longlst[k]; olati=latilst[k];
                 break;
              }
           dist_azimuth( latilst[i], longlst[i], olati, olong, &dist, &alpha);
           if( dist < 100 || dist > 1500 ) continue;
           sprintf(filename,"%s_amp_snr\0",sta_pair[j]);
           if((fin = fopen(filename, "r"))==NULL) {
//             printf("Amp file %s not found. Skipped\n",filename);
             continue;
           }
           for(k=0;;){
              if( fgets(buff, 300, fin) == NULL ) break;
              if((sscanf(buff,"%f %f %f %f %f", &period, &amppos, &snrpos, &ampneg, &snrneg))!=5) { 
                 cout<<"Wrong format! Skipped: "<<filename<<endl; 
                 break; 
              }
              if( period < perl ) continue;
              if( period > perh ) break;
              if( snrpos>5 || (snrpos>3 && snrneg>5) ) {
                 pertmp[k] = (period-per)*(period-per)+1e-10;
                 amptmp[k] = amppos;
                 k++;
              }
           }
           fclose(fin);
        }
        else if(strcmp(stalst[i],sta2[j])==0) {
           strcpy(osta,sta1[j]);
           for(k=0;k<nsta;k++)
              if(strcmp(osta,stalst[k])==0) {
                 olong=longlst[k]; olati=latilst[k];
                 break;
              }
           dist_azimuth( latilst[i], longlst[i], olati, olong, &dist, &alpha);
           if( dist < 100 || dist > 1500 ) continue;
           sprintf(filename,"%s_amp_snr\0",sta_pair[j]);
           if((fin = fopen(filename, "r"))==NULL) {
//             printf("Amp file %s not found. Skipped\n",filename);
             continue;
           }
           for(k=0;;){
              if( fgets(buff, 300, fin) == NULL ) break;
              if((sscanf(buff,"%f %f %f %f %f", &period, &amppos, &snrpos, &ampneg, &snrneg))!=5) { 
                 cout<<"Wrong format! Skipped: "<<filename<<endl;
                 break; }
              if( period < perl ) continue;
              if( period > perh ) break;
              if( snrneg>5 || (snrpos>5 && snrneg>3) ) {
                 pertmp[k] = (period-per)*(period-per)+1e-10;
                 amptmp[k] = ampneg;
                 k++;
              }
           }
           fclose(fin);
        }
        else continue;

        if( k<2 ) continue;
           for(iper=0;iper<2;iper++)
              for(jper=iper+1;jper<k;jper++) {
                 if(pertmp[iper]>pertmp[jper]) {
                    period=pertmp[iper];
                    pertmp[iper]=pertmp[jper];
                    pertmp[jper]=period;
                    amppos=amptmp[iper];
                    amptmp[iper]=amptmp[jper];
                    amptmp[jper]=amppos;
                 }
              }
        amppos=(amptmp[0]/pertmp[0]+amptmp[1]/pertmp[1])/(1./pertmp[0]+1./pertmp[1]);
        deg_l=(int)alpha;
        deg_h=deg_l+1;
        if(deg_h==360) deg_h=0;
        d_alpha=alpha-deg_l+0.5;
        amp_n[deg_l][0]+=amppos/pow(e_n,-1e-3*dist)/d_alpha;
        amp_n[deg_l][1]+=1./d_alpha;
        amp_n[deg_l][3]+=1;
        d_alpha=deg_h-alpha+0.5;
        amp_n[deg_h][0]+=amppos/pow(e_n,-1e-3*dist)/d_alpha;
        amp_n[deg_h][1]+=1./d_alpha;
        amp_n[deg_h][2]+=1;
/*        if(deg_l==204) cout<<"ampsta: "<<amppos<<"  amp_h: "<<amp_n[deg_h][0]<<"  weight: "<<amp_n[deg_h][1]<<"  num: "<<amp_n[deg_h][2]<<endl;
        if(deg_l==205) cout<<"ampsta: "<<amppos<<"  amp_l: "<<amp_n[deg_l][0]<<"  weight: "<<amp_n[deg_l][1]<<"  num: "<<amp_n[deg_l][3]<<endl; */
     }
     sprintf(outname,"%s_amp_azimuth",stalst[i]);
     if((fout = fopen(outname, "w"))==NULL) {
       printf("Cannot open file %s to write\n",outname);
       continue;
     }
     for(deg_l=0;deg_l<360;deg_l++)
        if(amp_n[deg_l][2]>0 && amp_n[deg_l][3]>0) fprintf(fout,"%3d   %.4g\n", deg_l, amp_n[deg_l][0]/amp_n[deg_l][1]);
     fclose(fout);
     cout<<"The "<<i<<"th station "<<stalst[i]<<" completed..."<<endl;
  }
  return 1;
}
