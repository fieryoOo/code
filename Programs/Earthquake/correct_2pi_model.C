#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;


double get_dist(double lat1,double lon1,double lat2,double lon2)
{
  double theta,pi,temp;
  double radius=6371;
  pi=4.0*atan(1.0);

  lat1=atan(0.993277*tan(lat1/180*pi))*180/pi;
  lat2=atan(0.993277*tan(lat2/180*pi))*180/pi;

  temp=sin((90-lat1)/180*pi)*cos(lon1/180*pi)*sin((90-lat2)/180*pi)*cos(lon2/180*pi)+sin((90-lat1)/180*pi)*sin(lon1/180*pi)*sin((90-lat2)/180*pi)*sin(lon2/180*pi)+cos((90-lat1)/180*pi)*cos((90-lat2)/180*pi);
  if(temp>1)
    {
      cout<<"warning cos(theta)>1 and correct to 1!!"<<temp<<endl;
      temp=1;
    }
  if(temp<-1)
    {
      cout<<"warning cos(theta)<-1 and correct to -1!!"<<temp<<endl;
      temp=-1;
    }
  theta=fabs(acos(temp));
  return theta*radius;
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



int main(int na, char *arg[])
{ 
  if(na!=3)
    {
      cout<<"usage:correct_2pi_model [in_file] [period]"<<endl;
      return 0;
    }
  FILE *ff,*f1,*f2;
  char path[300], mod_name[300], buff[300];
  double dis, distmp, perl, perh, per=atof(arg[2]);
  if((f1 = fopen(arg[1], "r"))==NULL) {
    printf("cannot open %s\n",arg[1]);
    exit(1);
  }
  int N=3000;
  char evnm[20], stnm[N][6];
  double lon[N],lat[N],time[N],time2[N],vel[N],amp[N],clon,clat,dtmp,temp,vell, velh;
  short YesNo[N], flag;
  int i,j,k,ii,jj,iev,ist;
  // cout<<"test"<<endl;
  fgets(buff, 300, f1);
  sscanf(buff,"%lf %lf %lf %lf %lf %s",&clon,&clat,&dtmp,&dtmp,&dtmp,evnm);
  sprintf(path,"/home/tianye/Model/PHASE_map/Event_Pre/%s/%s",evnm,evnm);
  for(i=0;i<N;i++)
    {
      if(fgets(buff,300,f1) == NULL) break;
      sscanf(buff,"%lf %lf %lf %lf %lf %s",&lon[i],&lat[i],&time[i],&vel[i],&amp[i],&stnm[i][0]);
      YesNo[i]=1;
    }
  ist=i;
  fclose(f1);

//  sprintf(buff,"%s_misfit",arg[1]);
//  ff=fopen(buff,"w");
//  sprintf(buff,"%s_vel_mod",arg[1]);
//  f2=fopen(buff,"w");
  for(i=0;i<ist;i++) {
     sprintf(mod_name,"%s_%s.dat",path,stnm[i]);
     if((f1=fopen(mod_name,"r"))==NULL) {
        cout<<"Cannot find mod file: "<<mod_name<<". Skipped!"<<endl;
        continue;
     }
     fgets(buff, 300, f1);
     perl=-1; flag=1;
     for(;;) {
        if(fgets(buff, 300, f1)==NULL) { flag=0; break; }
        sscanf(buff,"%lf %lf",&perh,&velh);
        if(fabs(perh-per)<=0.1) break;
        else if(perh>per && perl>0) {
           velh = (velh-vell)/(perh-perl)*(per-perl)+vell;
           break;
        }
        else if(perh>per) { flag=0; break; }
        perl=perh; vell=velh;
     }
     fclose(f1);
     if(flag==0) {
        cout<<"Period "<<per<<" is out of the modle range. Skipped!"<<endl;
        continue;
     }
     dis=vel[i]*time[i];
     //calc_dist(clat,clon,lat[i],lon[i],&distmp);
     //cout<<dis-distmp<<endl;
     temp=dis/velh;
     //fprintf(ff,"%lf %lf %lf\n",lon[i],lat[i],time[i]-temp);
     //time[i] -= 4.06;
     for(;;) {
        if(time[i]-temp>per*0.5) //0.25
           time[i]-=per;
        else if(time[i]-temp<-per*0.5) //0.75
           time[i]+=per;
        else break;
     }
     time2[i] = time[i];
     if(time2[i]-temp<-0.*per) time2[i]+=per*0.5; // -0.25
     //if(fabs(time2[i]-temp)>per/6.) { YesNo[i]=-1; std::cerr<<"threw!"<<std::endl; }
//cout<<perl<<" "<<per<<" "<<perh<<endl;
     vel[i]=dis/time2[i];
//     if(fabs(time[i]-temp)>per/4.) YesNo[i]=-1; 
//cout<<"misfit: "<<fabs(time[i]-temp)*2./per<<" pi"<<endl;}
//     fprintf(ff,"%lf %lf %lf\n",lon[i],lat[i],-(time[i]-temp)*2./per);
//     fprintf(f2,"%lf %lf %lf\n",lon[i],lat[i],velh);
  }
//  fclose(ff);
//  fclose(f2);

  sprintf(buff,"%s_v1",arg[1]);
  ff=fopen(buff,"w");
  fprintf(ff,"%lf %lf 0.0 NaN NaN %s\n",clon,clat, evnm);
  for(i=0;i<ist;i++) {
     if(YesNo[i]==1)
        fprintf(ff,"%lf %lf %lf %lf %lf %s\n",lon[i],lat[i],time[i],vel[i],amp[i], stnm[i]);
  }
  fclose(ff);

  sprintf(buff,"%s_v1_crctd",arg[1]);
  ff=fopen(buff,"w");
  fprintf(ff,"%lf %lf 0.0 %s\n",clon,clat, evnm);
  for(i=0;i<ist;i++) {
     if(YesNo[i]==1)
        fprintf(ff,"%lf %lf %lf %s\n",lon[i],lat[i],time2[i], stnm[i]);
  }
  fclose(ff);

  return 1;
}
