#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "/home/tianye/code/Programs/head/mysac.h"

using namespace std;

double get_dist2(double lat1, double lon1, double lat2, double lon2)
{
double pi;
pi = 4.0*atan(1.0);
double cva = 6378.137;
double cvb = 6356.7523142;
double f = 1/298.257223563;
double L = 0.00;
L = lon1-lon2;
//cout<<"L "<<L<<endl;
if (L > 180.000)  L =360.000000 - L;
if (L < -180.000) L =  360.000 - fabs(L);
L = fabs(L);
//cout<<"L: "<<L<<endl;
double U1 = 0;
U1 = atan((1-f)*tan(lat1/180*pi));
double U2 = 0;
U2 = atan((1-f)*tan(lat2/180*pi));
//cout<<"U1 "<<U1<<"U2 "<<U2<<" "<<sin(U1)<<" "<<sin(U2)<<" "<<cos(U1)<<" "<<cos(U2)<<endl;
double cv,cv1,cv2,cv3,cv4,cv5,cvC,numda1;
L = L*pi/180;
double numda = L;
numda1 = numda;
//cout<<"numda "<<numda<<"cos numda "<<cos(numda)<<endl;
do {
numda = numda1;
cv1 =  sqrt( (cos(U2)*sin(numda))*(cos(U2)*sin(numda))+ (cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda))*(cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda)) ); // cv1 sin(quan)
cv2 = sin(U1)*sin(U2)+ cos(U1)*cos(U2)*cos(numda);
//cout<<"cv1 cv2: "<<cv1<<" "<<cv2<<endl;
cv = atan2(cv1,cv2);
if(cv==0)cv+=0.00000000001;
cv3 = cos(U1)*cos(U2)*sin(numda)/sin(cv);
cv4 = 1 - cv3*cv3;
if(cv4==0)cv4+=0.00000000001;
cv5 = cos(cv) - 2*sin(U1)*sin(U2)/cv4;
cvC = f/16*cv4*(4 + f*(4 - 3*cv4));
numda1 = L + (1-cvC)*f*cv3*(cv + cvC*cv1*(cv5 + cvC*cv2*(-1 +2*cv5*cv5)));
//cout<<"numda1 "<<numda1<<endl;
} while (fabs(numda - numda1) > 0.0000000001);

double mius, cvA, cvB, deltacv,s;
mius = cv4*(cva*cva - cvb*cvb)/(cvb*cvb);
cvA = 1+mius/16384*(4096 + mius*(-768 + mius*(320 - 175*mius)));
//cout<<"cvA "<<cvA<<endl;
cvB = mius/1024*(256+ mius*(-128 + mius*(74 - 47*mius)));
deltacv = cvB*cv1*(cv5 +cvB/4*(cv2*(-1 + 2*cv5*cv5)-cvB/6*cv5*(-3+4*cv1*cv1)*(-3+4*cv5*cv5) ));
//cout<<"delatacv "<<deltacv<<"cvb "<<cvb<<"cv "<<cv<<endl;
s = cvb * cvA *(cv - deltacv);
return s;
}


int main (int argc, char *argv[])
{

  FILE *f,*flst,*ff,*fout;
  char outname[300],buff[300];
  char sta[500][8],filelist[5000][100];
  double dislong,dislati,dist,longout;
  double stalong[500],stalati[500],data[500],deg[500];
  int i,ii,ista,ilst,ifile,nfile;

  int flag,ibin,sta_num[8];
  double bin_avg[8],bin_weight[8],avg,weight;
  double PI=4*atan2(1.,1.);
  if(argc != 3)
    {
      printf("Usage: bin_avg_center_sta.C [input_file_lst each_centersta_1st_line] [sta.lst]\n");
      exit(-1);
    }

  if((f = fopen(argv[2],"r"))==NULL){
     fprintf(stderr,"cannot open sta_lst file %s\n",argv[2]);
     return 0;
    }

  char stalst[5000][8];
  double longlst[5000],latilst[5000];
  for (i=0;;i++) {
     if (fscanf(f,"%s %lf %lf", &stalst[i],&longlst[i],&latilst[i]) != 3) break;
    }
  ilst=i;
  fclose(f);


  if((flst = fopen(argv[1],"r"))==NULL){
     fprintf(stderr,"cannot open data file list %s\n",argv[1]);
     return 0;
    }

  for (i=0;;i++) {
    if (fscanf(flst,"%s", &filelist[i]) != 1) break;
    }
  nfile = i;
  fclose(flst);

//  sprintf(outname,"%s.txt",argv[1]);
//  if((fout = fopen(outname, "wb"))==NULL) {
//     fprintf(stderr,"Could not open %s to write\n", outname);
//     return 0;
//    }

 for (ifile=0;ifile<nfile;ifile++) {
//printf("%d: %s\n",ifile,filelist[ifile]);
  if((ff = fopen(filelist[ifile],"r"))==NULL){    
     fprintf(stderr,"cannot open data file %s\n",filelist[ifile]);
     continue;
    }

  data[0]=NAN;
  for (i=0;;i++) {
     if (i==0) fscanf(ff,"%s %lf %lf", &sta[i], &stalong[i], &stalati[i]);
     else if (fscanf(ff,"%s %lf", &sta[i], &data[i]) != 2)break;
//printf("%s %f %f\n",sta[i],stalong[i],stalati[i]);
    }
  ista=i;
  fclose(ff);

  sprintf(outname,"%s.txt",filelist[ifile]);
  if((fout = fopen(outname, "wb"))==NULL) {
       fprintf(stderr,"Could not open %s to write\n", outname);
     return 0;
    }

//printf("%d %d\n",ilst,ista);
  deg[0]=0;
  for (i=0;i<ilst;i++) { 
     for (ii=1;ii<ista;ii++){
//printf("sta: %s  lst_sta: %s\n",sta[ii],stalst[i]);
        if (strcmp(sta[ii],stalst[i])==0){
           stalong[ii]=longlst[i];
           stalati[ii]=latilst[i];
           dislong=get_dist2(stalati[0],stalong[0],stalati[0],stalong[ii]);
//printf("long1: %f  long2: %f  dislong: %f\n",stalong[0],stalong[ii],dislong);
           if(stalong[ii]<stalong[0])dislong=-dislong;
           dislati=get_dist2(stalati[0],stalong[ii],stalati[ii],stalong[ii]);
           if(stalati[ii]<stalati[0])dislati=-dislati;
           deg[ii] = atan2(dislati,dislong) * 180 / PI + 180;
//printf("dislati: %f  dislong: %f  deg: %f\n",dislati,dislong,deg[ii]);
          }
       }
    }

  flag=0;
//  avg=0;
//  weight=0;

  for (i=0;i<8;i++){sta_num[i]=0;}
  for (ii=1;ii<ista;ii++){
     ibin=(int)(deg[ii]/45.0);
//printf("deg: %f\n",deg[ii]);
     if(deg[ii]==360)ibin=ibin-1 ;
//     dist=get_dist2(stalati[0],stalong[0],stalati[ii],stalong[ii]);
       sta_num[ibin]=sta_num[ibin]+1;
//     bin_weight[ibin]+=1/dist;
//     bin_avg[ibin]+=data[ii]/dist;
//printf("deg: %lf  ibin: %d\n",deg[ii],ibin);
    }
  for (i=0;i<8;i++){
//printf("sta_num: %d\n",sta_num[i]);
     if(sta_num[i]==0)continue;
     flag++;
//     bin_avg[i]=bin_avg[i]/sta_num[i];
//     bin_weight[i]=bin_weight[i]/sta_num[i];
//printf("bin_avg: %f  sta_num: %d\n",bin_avg[i],sta_num[i]);
//     avg+=bin_avg[i];
//     weight+=bin_weight[i];
    }
  if(flag<=7){
     printf("Station %s: Got coverage in only %d out of 8 directions, skipped!\n",sta[0],flag);
     sprintf(buff,"rm -f %s\n",outname);
     system(buff);
    }
  else {
     printf("Station %s: Coverage in %d out of 8 directions.\n",sta[0],flag);
//     longout=stalong[i];
//     if(longout<0)longout+=360;
     for (ii=0;ii<ista;ii++){
       longout=stalong[ii];
       if(longout<0)longout+=360;
       fprintf(fout,"%f %f %f\n",longout,stalati[ii],data[ii]);
       }
    }
  fclose(fout);
 }
 return 0;
}
