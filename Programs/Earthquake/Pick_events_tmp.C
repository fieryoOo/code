#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;
#define NEV 20000

int calc_jday ( int y, int m, int d )
        /*------------------------------------------------------------------------*/
{
        int jd = 0;
        int i;

        for ( i = 1; i < m; i++ )
        {
                if ( (i==1) || (i==3) || (i==5) || (i==7) || (i==8) || (i==10) ) jd += 31;
                else if (i==2)
                {
                        if ( (y%400 == 0) || (y%4 == 0 && y%100 != 0 ) ) jd += 29;
                        else jd += 28;
                }
                else jd += 30;
        }

        return jd + d;
}


double abs_time ( int yy, int jday, int hh, int mm, int ss, int ms )
        /*--------------------------------------------------------------------------
          computes time in s relative to 1900
          --------------------------------------------------------------------------*/
{
        int nyday = 0, i;
        double abssec;
        for ( i = 1901; i < yy; i++ )
        {
                if ( (i%400==0) || (i%4==0 && i%100!=0) ) nyday += 366;
                else nyday += 365;
        }
        abssec =  24.*3600.*(nyday+jday) + 3600.*hh + 60.*mm + 1.*ss + 0.001*ms;
        return abssec;
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


main(int na, char *arg[])
{
  if (na!=3)
    {
      cout<<"usage: File_Combiner [ev_infile (ev lon lat ctlg XXX)] [ev_outfile]"<<endl;
      return 0;
    }

  FILE *fl;
  int i, j, ii, iflag, nev, nrms, nrmc, flag[NEV];
  int yer, mth, day, hor, min, sec, jday;
  double lon[NEV], lat[NEV], abtime[NEV], Adist[NEV], ftmp, dist;
  char info[NEV][300], ctlg[NEV][10], buff[300];

  if((fl=fopen(arg[1],"r"))==NULL) {
     cout<<"Cannot open file: "<<arg[1]<<endl;
     exit(0);
  }
  for(i=0;;i++) {
     if((fgets(info[i], 300, fl))==NULL) break;
     strtok(info[i],"\n");
     sscanf(info[i],"%4d%2d%2d%2d%2d%2d %lf %lf %s %lf %lf %lf",&yer,&mth,&day,&hor,&min,&sec,&lon[i],&lat[i],&ctlg[i],&ftmp,&ftmp,&Adist[i]);
     //cout<<yer<<" "<<mth<<" "<<day<<" "<<hor<<" "<<min<<" "<<sec<<" "<<lon[i]<<" "<<lat[i]<<endl;
     jday= calc_jday(yer,mth,day);
     abtime[i] = abs_time(yer, jday, hor, min, sec, 0 );
  }
  nev=i;
  fclose(fl);

  for(i=0;i<nev;i++){
     ftmp=abtime[i]; ii=i;
     for(j=i+1;j<nev;j++)
        if(ftmp>abtime[j]) {
           ftmp=abtime[j];
           ii=j;
        }
     if(ii==i) continue;
     ftmp=abtime[i];
     abtime[i]=abtime[ii];
     abtime[ii]=ftmp;
     ftmp=lon[i];
     lon[i]=lon[ii];
     lon[ii]=ftmp;
     ftmp=lat[i];
     lat[i]=lat[ii];
     lat[ii]=ftmp;
     ftmp=Adist[i];
     Adist[i]=Adist[ii];
     Adist[ii]=ftmp;
     sprintf(buff,"%s",info[i]);
     sprintf(info[i],"%s",info[ii]);
     sprintf(info[ii],"%s",buff);
     sprintf(buff,"%s",ctlg[i]);
     sprintf(ctlg[i],"%s",ctlg[ii]);
     sprintf(ctlg[ii],"%s",buff);
  }

  nrms=0; nrmc=0;
  for(i=0;i<nev;i++) flag[i]=1;
  for(i=0;i<nev;i++) {
     ii=0; iflag=flag[i];
     for(j=i+1;j<nev;j++) {
        if(abtime[j]-abtime[i]>500.) break;
        //calc_dist(lat[i], lon[i], lat[j], lon[j], &dist);
        if(fabs(Adist[i]-Adist[j])<2000.) {
           cout<<abtime[j]-abtime[i]<<":  "<<info[i]<<"  -  "<<info[j]<<endl;
           flag[i]=0;
           if(strcmp(ctlg[i],ctlg[j])==0) { ii+=1;  nrmc+=flag[j]; flag[j]=0; }
        }
     }
     if(ii>0) nrmc+=iflag;
     else nrms+=iflag-flag[i];
  }
  cout<<nrms+nrmc<<" events removed: "<<nrms<<" repeated and "<<nrmc<<" merged."<<endl;

  fl=fopen(arg[2],"w");
  for(i=0;i<nev;i++)
     if(flag[i]==1) fprintf(fl,"%s\n",info[i]);
  fclose(fl);

  return 1;

}
