#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

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

  FILE *f,*fout;
  char outname[300],buff[300];
  double long0,lati0,dislong,dislati,avg,weight,weightsum,delta;
  double stalong[5000],stalati[5000],data[5000],deg[5000],dist[5000],azidata[1000],azideg[1000],azidist[1000];
  int i,ii,ista,iazi,imax,imin;

  double PI=4*atan2(1.,1.);

  if(argc != 4)
    {
      printf("Usage: amp_dist_angle.C [center_sta_long] [center_sta_lati] [input_data_file]\n");
      exit(-1);
    }

  long0=atof(argv[1]);
  lati0=atof(argv[2]);

  if((f = fopen(argv[3],"r"))==NULL){
     fprintf(stderr,"cannot open data file %s\n",argv[3]);
     return 0;
    }
  for (i=0;;i++) {
     if (fscanf(f,"%lf %lf %lf", &stalong[i],&stalati[i],&data[i]) != 3) break;
    }
  ista=i;
  fclose(f);

  sprintf(outname,"%s_azi",argv[3]);
  if((fout = fopen(outname, "wb"))==NULL) {
       fprintf(stderr,"Could not open %s to write\n", outname);
     return 0;
    }

  for (i=0;i<ista;i++) {
     dislong=get_dist2(lati0,long0,lati0,stalong[i]);
     if(stalong[i]<long0)dislong=-dislong;
     dislati=get_dist2(lati0,stalong[i],stalati[i],stalong[i]);
     if(stalati[i]<lati0)dislati=-dislati;
     deg[i] = atan2(dislati,dislong) * 180 / PI + 180;
     dist[i]=get_dist2(lati0,long0,stalati[i],stalong[i]);
     if(dist[i]>100 && dist[i]<500){
        fprintf(fout,"%f %f %f\n",deg[i],log(data[i]),dist[i]);
        if(deg[i]<10)fprintf(fout,"%f %f %f\n",deg[i]+360,log(data[i]),dist[i]);
        if(deg[i]>350)fprintf(fout,"%f %f %f\n",deg[i]-360,log(data[i]),dist[i]);
       }
     }
  fclose(fout);

  sprintf(buff,"sort -g %s > temp_azi",outname);
  system(buff);
  sprintf(buff,"mv temp_azi %s",outname);
  system(buff);

  if((f = fopen(outname,"r"))==NULL){
     fprintf(stderr,"cannot open azi file %s\n",outname);
     return 0;
    }
  delta=1.;
  for (i=0;;i++) {
     if (fscanf(f,"%lf %lf %lf", &azideg[i],&azidata[i],&azidist[i]) != 3) break;
    }
  iazi=i;
  fclose(f);

  sprintf(outname,"%s_azi_factor",argv[3]);
  if((fout = fopen(outname, "wb"))==NULL) {
       fprintf(stderr,"Could not open %s to write\n", outname);
     return 0;
    }

  for (i=0;i<360;i++){
      weightsum=0;avg=0;
      for (ii=0;ii<iazi;ii++){
          if(azideg[ii]>i-4 && azideg[ii]<i+4){
             weight=exp(-(i-azideg[ii])*(i-azideg[ii])/2/delta/delta);
             avg+=azidata[ii]*weight;
             weightsum+=weight;
            }
         }
      avg=avg/weightsum;
      fprintf(fout,"%d %f\n",i,avg);
     }
  fclose(fout);
  return 0;
}
