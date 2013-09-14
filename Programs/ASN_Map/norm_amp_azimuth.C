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
  
  FILE *fin, *fout;
  char filename[300], outname[300], buff[300];
  char stalst[5000][10], *tmp, *staname;
  int flag[3600], flagl[3600]={0}, flagh[3600]={0};
  int i, j, k, ii, namp, iper, jper, nfile, nsta, nost;
  float olong2[3000], olati2[3000];
  float olong[5000], olati[5000], ampst[5000], longlst[5000], latilst[5000];
  float per, perl, perh, pertmp[30], amptmp, fstep=-0.08;
  float amppos, ampneg, clong, clati, longtmp, latitmp;
  double dist, alpha[5000], alphatmp, powe, amp_n[3600][100]={0}, weight;
  double amp_tmp[2000],avg,std;
  double e_n=2.7182818284590451;
  if(argc != 5) 
    {
      printf("Usage: norm_amp_azimuth [input_map] [center long] [center lati] [min dist]\n"); 
      exit(-1);
    }
  
  for(i=0;i<3600;i++) amp_n[i][0]=0;
  clong=atof(argv[2]);
  if(clong<0) clong+=360;
  else if(clong>360) clong-=360;
  clati=atof(argv[3]);
  if(clati>90 || clati<-90) {cout<<"Wrong latitude!"<<endl; exit(0);}
  if((fin = fopen(argv[1], "r"))==NULL) {
    printf("Cannot open input file %s\n",argv[1]);
    exit(1);
  }
  for(i=0;;i++) {
    if((fscanf(fin,"%f %f %f", &olong[i], &olati[i], &ampst[i]))!=3) break;
    if(olong[i]<0) olong[i]+=360;
    else if(olong[i]>360) olong[i]-=360;
  }
  nsta = i;
  fclose(fin);

  fout = fopen("temp.txt", "w");
  nost=0;
  for(i=0;i<nsta;i++) {
     calc_dist( clati, clong, olati[i], olong[i], &dist );
     if( dist < atof(argv[4])-20 ) { ampst[i]=-1;continue; }
     if( dist > atof(argv[4])+500 ) continue;
     olong2[nost]=olong[i];
     olati2[nost]=olati[i];
     calc_azimuth( clati, clong, olati2[nost], olong2[nost], &alpha[nost]);
     fprintf(fout,"%f %f %f\n",olong[i],olati[i],ampst[i]/pow(e_n,-1e-3*dist));
     nost++;
  }
  fclose(fout);
  if(nost==0) exit(0);

  for(k=0;k<20;k++)
     for(j=0;j<nost;j++) {
        if(alpha[j]*10<k+20 && alpha[j]*10>=k) flagh[k]=1;
        else if(alpha[j]*10<=k || alpha[j]*10>k+3580) flagl[k]=1;
     }
  for(k=20;k<3581;k++)
     for(j=0;j<nost;j++) {
        if(alpha[j]*10<k+20 && alpha[j]*10>=k) flagh[k]=1;
        else if(alpha[j]*10<=k && alpha[j]*10>k-20) flagl[k]=1;
     }
  for(k=3581;k<3600;k++)
     for(j=0;j<nost;j++) {
        if(alpha[j]*10<k-3580 || alpha[j]*10>=k) flagh[k]=1;
        else if(alpha[j]*10<=k && alpha[j]*10>k-20) flagl[k]=1;
     }
  for(k=0;k<3600;k++) flag[k]=flagh[k]*flagl[k];

  sprintf(buff,"minmax temp.txt | awk '{print $5,$6}' | sed s/'<'/''/g | sed s/'>'/''/g | sed s/'\\/'/' '/g | awk '{printf \"-R%%.0f/%%.0f/%%.0f/%%.0f\\n\",$1-3,$2+3,$3-3,$4+3}' > region_temp");
  system(buff);
  sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel_positive temp.txt region_temp 0.1");
  system(buff);

  if((fin = fopen("temp.txt.HD", "r"))==NULL) {
    printf("Cannot open file temp.txt.HD to read\n");
    exit(0);
  }
  for(j=0;;j++) {
     if( fgets(buff, 300, fin) == NULL ) break;
     if((sscanf(buff,"%f%f%f", &longtmp, &latitmp, &amppos))!=3) {
        cout<<"Wrong format! Stopped: "<<filename<<endl;
        break;
     }
     calc_dist( clati, clong, latitmp, longtmp, &dist );
     if(dist<atof(argv[4])+50 || dist>atof(argv[4])+350) continue;
     calc_azimuth( clati, clong, latitmp, longtmp, &alphatmp );
     alphatmp=(int)(alphatmp*10+0.5);
     if(alphatmp==3600) alphatmp=0;
     if(flag[(int)alphatmp]==0) continue;
     amp_n[(int)alphatmp][0]+=1;
     amp_n[(int)alphatmp][(int)amp_n[(int)alphatmp][0]]=amppos;
  }
  fclose(fin);
/*
  for(k=0;k<3600;k++) {
     if(flag[k]==0 || amp_n[k][0]==0) {amp_n[k][1]=-1; cout<<"0"<<endl;}
     else {amp_n[k][1]=amp_n[k][1]/amp_n[k][0]; cout<<amp_n[k][0]<<endl;}
//     cout<<"amp: "<<amp_n[k][0]<<"  snum: "<<amp_n[k][1]<<endl;
  }
*/
//  cout<<amp_n[0][0]<<endl;
  int jj = 0;
  int cvk = 0;
  for(k=0;k<3600;k++){
     weight=0; namp=0; avg=0;
//     cout<<amp_n[0][0]<<endl;
     cvk = 0;

     for(j=k-30;j<=k+30;j++){
 //       cout<<k<<" "<<j<<endl;
        jj = j;
	if (j<0) jj = j+3600;
	if (j>=3600) jj = j - 3600;
 // 	cout<<"check jj "<<jj<<"amp_n "<<amp_n[jj][0]+1<<endl;
//        cout<<"cv "<<amp_n[jj][0]<<" "<<jj<<endl;

        for(ii=1;ii<amp_n[jj][0]+1;ii++) {
//	   cout<<"ii "<<ii<<" namp "<<namp<<" "<<namp+ii-1<<endl;
           amp_tmp[cvk]=amp_n[jj][ii];
           avg+=amp_n[jj][ii]*pow(e_n,-(j-k)*(j-k)/200.);
           cvk ++;
           }

        weight+=amp_n[jj][0]*pow(e_n,-(j-k)*(j-k)/200.);
// amp_tmp[61],avg,std;
     }
     cout<<cvk<<endl;
     amp_n[k][99]=avg/(weight+1e-20);
     std=0;
     for(ii=0;ii<cvk;ii++){
        std += pow((amp_tmp[ii]-amp_n[k][99]),2);
     }
     if (cvk<2)
        amp_n[k][98] = 50000;
     else 
        amp_n[k][98] = sqrt(std/(cvk-1));
     cout<<amp_n[k][98]<<endl;
  }
/*  
  fout=fopen("temp_azimuth","w");
  for(k=0;k<3600;k++) {
     if(amp_n[k][1]==0)continue;
     fprintf(fout,"%f  %g  %g\n",k/10.,amp_n[k][1],amp_n[k][0]);
  }
  fclose(fout);
*/

  for(i=0;i<nsta;i++){
     if(ampst[i]==-1) continue;
     calc_azimuth( clati, clong, olati[i], olong[i], &alphatmp);
     alphatmp=(int)(alphatmp*10+0.5);
     if(alphatmp==3600) alphatmp=0;
     if(amp_n[(int)alphatmp][1]==0){ampst[i]=-1;continue;}
     ampst[i]=ampst[i]/amp_n[(int)alphatmp][1];
  }

  sprintf(filename,"%s_azi\0",argv[1]);
  fout=fopen(filename,"w");
  for(k=0;k<3600;k++) {
     if(amp_n[k][1]==0)continue;
     fprintf(fout,"%f  %g %g\n",k/10.,amp_n[k][99],amp_n[k][98]);
  }
  fclose(fout);

  sprintf(filename,"%s_azi_norm",argv[1]);
  if((fout=fopen(filename,"w"))==NULL) {
     cout<<"Can't open file "<<filename<<" to write!";
     exit(0);
  }
  for(i=0;i<nsta;i++)if(ampst[i]>=0)
     fprintf(fout,"%.3f  %.3f  %.5g\n",olong[i],olati[i],ampst[i]);
  fclose(fout);
}
