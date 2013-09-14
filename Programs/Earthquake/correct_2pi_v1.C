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
  if(na!=4)
    {
      cout<<"usage:correct_2pi event_file_name period least_sta_#"<<endl;
      return 0;
    }
  FILE *ff,*f1,*f2;
  char file_name[400];
  double dis,dis_min;
  double per;
  per=atof(arg[2]);
  if((f1 = fopen(arg[1], "r"))==NULL) {
    printf("cannot open %s\n",arg[1]);
    exit(1);
  }
  cout<<"Correcting 2pi for file "<<arg[1]<<"..."<<endl;
  int N=3000;
  char sta[N][10], evbuff[300];
  double lon[N],lat[N],time[N],vel[N],amp[N],clon,clat,dtmp;
  int YesNo[N];
  int i,j,k,ii,jj,iev,ist;
  // cout<<"test"<<endl;
  fgets(evbuff, 300, f1);
  sscanf(evbuff,"%lf %lf %lf %lf %lf",&clon,&clat,&dtmp,&dtmp,&dtmp,&sta);
  for(i=0;i<N;i++)
    {
      if(fscanf(f1,"%lf %lf %lf %lf %lf %s",&lon[i],&lat[i],&time[i],&vel[i],&amp[i],&sta[i][0])==EOF) 
	break;
      YesNo[i]=0;
    }
  ist=i;

  //cout<<ist<<endl;
  fclose(f1);
  if(ist<atof(arg[3]))
    {
      printf("Only %d stations\n",ist);
      exit(1);
    }
  int order[ist];
  double fact[ist];
  for(i=0;i<ist;i++)
    {
      //      fact[i]=lat[i]-lon[i];
      fact[i]=(lat[i]-clat)*(lat[i]-clat)+(lon[i]-clon)*(lon[i]-clon);
    }
  order[0]=0;
  for(i=1;i<ist;i++)
    {
      for(j=0;j<i;j++)
	{
	  //	  if(fact[i]>fact[order[j]])
	  if(fact[i]<fact[order[j]])
	  break;
	}
      for(k=i;k>j;k--)
	{
	  order[k]=order[k-1];
	}
      order[j]=i;      
    }
  int mark,temp, nmsft;
  char buff[300];
  sprintf(buff,"%s_v1",arg[1]);
  //ff=fopen(buff,"w");
  // fprintf(ff,"%lf %lf %lf %lf %lf\n",lon[order[0]],lat[order[0]],time[order[0]],vel[order[0]],amp[order[0]]);
  int ino=0;
  //  for(i=0;i<ist;i++)
  //cout<<i<<" "<<lon[order[i]]<<" "<<lat[order[i]]<<endl;

/*
  calc_dist(lat[order[0]],lon[order[0]],clat,clon,&dis);
  temp=dis/vel[order[0]];
cout<<temp<<" "<<time[order[0]]<<"!!\n\n"<<endl;
  for(;;) {
     if(time[order[0]]-temp>per/2)
        time[order[0]]-=per;
     else if(temp-time[order[0]]>per/2)
        time[order[0]]+=per;
     else break;
  }
*/

  for(i=1;i<ist;i++)
    {
      dis_min=999999999;
      for(j=0;j<i;j++)
	{
	  //if(j==137)
	    //	    fprintf(stderr,"%d %d %lf %lf\n",i,j,lon[order[j]],lat[order[j]]);
	  //	  cout<<lat[order[i]]<<" "<<lon[order[i]]<<" "<<lat[order[j]]<<" "<<lon[order[j]]<<endl;
	  //dis=get_dist(lat[order[i]],lon[order[i]],lat[order[j]],lon[order[j]]);
          calc_dist(lat[order[i]],lon[order[i]],lat[order[j]],lon[order[j]],&dis);
	  // cout<<i<<" "<<j<<" "<<dis<<endl;
	  if(dis<dis_min)
	    {
	      dis_min=dis;
	      mark=j;
	      // cout<<mark<<endl;
	    }
	  //	  if(i==326)
	  // fprintf(stderr,"GILL %d %d\n",i,j);
	}
//cout<<dis_min<<endl;
      //      if(i==326)
      //fprintf(stderr,"%d %d %lf %lf %lf %lf\n",i,mark);
      //      cout<<">"<<endl;
      //cout<<lon[order[i]]<<" "<<lat[order[i]]<<endl;
      //cout<<lon[order[mark]]<<" "<<lat[order[mark]]<<endl;
      if(fabs(amp[order[i]]-amp[order[mark]])>10*amp[order[i]]||fabs(amp[order[i]]-amp[order[mark]])>10*amp[order[mark]])
	{
	  //	  if(i==137)
	  // fprintf(stderr,"%d %d\n",i,mark);
	  //cout<<dis<<" "<<amp[order[i]]<<" "<<amp[order[mark]]<<endl;
	  YesNo[order[i]]=999;
	  ino++;
	  lon[order[i]]=lon[order[mark]];
	  lat[order[i]]=lat[order[mark]];
	  time[order[i]]=time[order[mark]];
	  vel[order[i]]=vel[order[mark]];
	  amp[order[i]]=amp[order[mark]];
	  continue;
	}
      dis=vel[order[i]]*time[order[i]];
      temp=(int)(dis/vel[order[mark]]+0.5);
      for(;;)
	{
	  if(time[order[i]]-temp>per*0.5)
	    time[order[i]]-=per;
	  else
	    if(temp-time[order[i]]>per*0.5)
	      time[order[i]]+=per;
	    else
	      break;
	}
      if(fabs(time[order[i]]-temp)>6)
      //if(fabs(time[order[i]]-temp)>6 && fabs(fabs(time[order[i]]-temp)-per/2.)>6)
	{
	  //	  if(i==137)
	  //fprintf(stderr,"%d %d %lf %lf\n",i,mark,lon[order[mark]],lat[order[mark]]);
	  //cout<<"misft > 6 sec"<<endl;
	  YesNo[order[i]]=999;
	  ino++;
	  lon[order[i]]=lon[order[mark]];
	  lat[order[i]]=lat[order[mark]];
	  time[order[i]]=time[order[mark]];
	  vel[order[i]]=vel[order[mark]];
	  amp[order[i]]=amp[order[mark]];
	  continue;
	}
      vel[order[i]]=dis/time[order[i]];
      //fprintf(ff,"%lf %lf %lf %lf %lf\n",lon[order[i]],lat[order[i]],time[order[i]],vel[order[i]],amp[order[i]]);
      //      fprintf(ff,"%lf %lf %lf %lf %lf %lf %lf %lf\n",lon[order[i]],lat[order[i]],time[order[i]],vel[order[i]],amp[order[i]],lon[order[mark]],lat[order[mark]],vel[order[mark]]);
      //cout<<lon[order[i]]<<" "<<lat[order[i]]<<" "<<fact[order[i]]<<endl;
    }
  //  cout<<ino<<" "<<ist<<" "<<ino/ist-0.9<<endl;
  if(ino*1.0/ist>0.3||ist-ino<atof(arg[3]))
    {
      printf("No result: %d out of %d stations removed.\n",ino,ist);
      exit(1);
    }
  ff=fopen(buff,"w");
 
  fprintf(ff,"%s",evbuff); 
  for(i=0;i<ist;i++)
    {
      if(YesNo[order[i]]==0)
	fprintf(ff,"%lf %lf %lf %lf %lf %s\n",lon[order[i]],lat[order[i]],time[order[i]],vel[order[i]],amp[order[i]],sta[order[i]]);
    }
  
  fclose(ff);
  cout<<"Result OK."<<endl;
  return 0;
}
