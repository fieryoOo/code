#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
using namespace std;

#define NSTA 2000

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


int gau_fit_line ( double *datax, double *datay, int idata, double slpmn, double slpmx, double alpha, double *itcptout )
{
  FILE *f;
  int i,j,ii;
  double weit, weight;
  double slp,slpstp,itcpt,itcptmx,itcptmn,itcptstp,slpfit,itcptfit,slpp[100],itcptt[100],rsd[100],dataxx,datayy,datamin,rsdmn,rsdfit;

  for(i=0;i<idata;i++) {
     datamin=datax[i]; ii=i;
     for(j=i+1;j<idata;j++)
        if(datamin>datax[j]){
           datamin=datax[j];
           ii=j;
        }
     if(ii==i) continue;
     dataxx=datax[i];
     datayy=datay[i];
     datax[i]=datax[ii];
     datay[i]=datay[ii];
     datax[ii]=dataxx;
     datay[ii]=datayy;
  }

  itcptmx=-999999999.;itcptmn=999999999.;
  for(i=0;i<idata;i++){
     itcpt=datay[i]-datax[i]*slpmx;
     if(itcptmn>itcpt)itcptmn=itcpt;
     itcpt=datay[i]-datax[i]*slpmn;
     if(itcptmx<itcpt)itcptmx=itcpt;
    }


  for(;;){
//    printf("slpmin: %f  slpmax: %f\n",slpmn,slpmx);
    slpstp=(slpmx-slpmn)/9;
    itcptstp=(itcptmx-itcptmn)/9;

    rsdmn=999999999.;
    for(i=0;i<10;i++){
       slp=slpmn+i*slpstp;
       for(j=0;j<10;j++){
          itcpt=itcptmn+j*itcptstp;
          rsd[i*10+j]=0; weight=0;
          slpp[i*10+j]=slp;
          itcptt[i*10+j]=itcpt;
          for(ii=0;ii<idata;ii++){
             weit=exp(-alpha*pow(datax[ii],2)); //weight+=weit;
             rsd[i*10+j] += pow((datay[ii]-itcpt-datax[ii]*slp),2)*weit;
             //rsd[i*10+j]+=(slp*datax[ii]+itcpt-datay[ii])*(slp*datax[ii]+itcpt-datay[ii]);
            }
          //rsd[i*10+j]=sqrt(rsd[i*10+j]/(idata-1));
          rsd[i*10+j] = sqrt(rsd[i*10+j]);///(weight-1);
         }
      }
    for(i=0;i<5;i++)
       for(j=i;j<100;j++)
          if(rsd[i]>rsd[j]){
             rsdfit=rsd[i];
             slpfit=slpp[i];
             itcptfit=itcptt[i];
             rsd[i]=rsd[j];
             slpp[i]=slpp[j];
             itcptt[i]=itcptt[j];
             rsd[j]=rsdfit;
             slpp[j]=slpfit;
             itcptt[j]=itcptfit;
            }

    slpmn=999999999;slpmx=-999999999;
    itcptmn=999999999;itcptmx=-999999999;
    for(i=0;i<5;i++){
//printf("residue: %f\n",rsd[i]);
       if(slpmn>slpp[i])slpmn=slpp[i];
       if(slpmx<slpp[i])slpmx=slpp[i];
       if(itcptmn>itcptt[i])itcptmn=itcptt[i];
       if(itcptmx<itcptt[i])itcptmx=itcptt[i];
      }
    if((slpmn-slpmx)*(slpmn-slpmx)<slpmn*slpmn/100000)break;
    //cout<<slpmx<<" "<<slpmn<<endl;
   }
//  printf("slope: %f  intercept: %f  residue: %f\n",slpp[0],itcptt[0],rsd[0]);
//  printf("%f %f %f\n",slpp[0],itcptt[0],rsd[0]);
//  *slpout=slpp[0]; *rsdout=rsd[0]; 
  *itcptout=itcptt[0];
  return 1;
}

int Gauss_Smoothing(char *fname, double hdis)
{
   FILE *f1;
   int i,j,nsta,ndata;
   char buff[300];
   double lon[NSTA], lat[NSTA], dat[NSTA];
   double dis, itcpt, dist[NSTA/4], data[NSTA/4];
   double alpha=0.5/hdis/hdis;
   if((f1=fopen(fname,"r")) == NULL) {
      cout<<"Cannot open file "<<fname<<endl;
      return 0;
   }
   for(nsta=0;;nsta++){
      if(fgets(buff, 300, f1) == NULL) break;
      sscanf(buff,"%lf %lf %lf", &lon[nsta], &lat[nsta], &dat[nsta]);
   }
   fclose(f1);

   double dsmd[nsta], mean, slpmn, slpmx;

   for(i=0;i<nsta;i++){
      ndata=0; mean=0;
      for(j=0;j<nsta;j++){
         calc_dist(lat[i], lon[i], lat[j], lon[j], &dis);
         if(dis>2*hdis) continue;
         dist[ndata]=dis; data[ndata]=dat[j];
         mean += dat[j];
         ndata++;
      }
      mean = fabs(mean/ndata);
      slpmx = 0.05*mean; slpmn = -slpmx;
      //f1=fopen("temp","w");
      //fprintf(f1,"%lf %lf\n",slpmn, slpmx);
      //for(j=0;j<ndata;j++) fprintf(f1,"%lf %lf\n",dist[j], data[j]);
      //fclose(f1);
      if(ndata<3) itcpt=dat[i];
      else gau_fit_line(&dist[0], &data[0], ndata, slpmn, slpmx, alpha, &itcpt);
      dsmd[i]=itcpt;
   }

   sprintf(buff,"%s_smd\0",fname);
   f1=fopen(buff,"w");
   for(i=0;i<nsta;i++)
      fprintf(f1,"%lf %lf %lf\n", lon[i], lat[i], dsmd[i]);
   fclose(f1);

   return 1;
}

int main(int argc, char *argv[])
{
   if(argc!=3){
      cout<<"Usage: Gauss_Smoothing [input_file] [half_dist]"<<endl;
      return 0;
   }

   double hdis=atof(argv[2]);
   Gauss_Smoothing(argv[1], hdis);

   return 1;
}
