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


int main(int na, char *arg[])
{
  if(na!=6)
    {
      cout<<"usage:correct_tr_t_curvature event_name period least_sta_# greatest_curvature amp_curvature"<<endl;
      return 0;
    }
  FILE *ff,*fin,*fin_amp,*fin2,*fin3,*fout,*file1;
  int i,j;
  int npts_x,npts_y;
  char buff1[300],sta1[10],path[300],stnm[10], evbuff[100];
  double lat,lon,lat2,lon2,t_lat,t_lon,radius,pi,sta1_lon,sta1_lat;
  int t_i,t_j;
  int marker_nn,marker_EN[2][2],marker_E,marker_N;
  double period,dist;
  period=atof(arg[2]);

  radius=6371.1391285;
  pi=4.0*atan(1.0);
  double dx,dy,x0,y0,x1,y1,temp,temp1,temp2,temp3,temp4,lat_temp;
//  npts_x=63;
//  npts_y=61;

  dx=0.1;//degree
  dy=0.1;//degree

 x0=999;
 x1=-999;
 y0=999;
 y1=-999;

  cout<<"Correcting curvature for event "<<arg[1]<<endl;
  sprintf(buff1,"%s.ph.txt_v1.HD",arg[1]);
  if((fin=fopen(buff1,"r"))==NULL)
    {
      cout<<buff1<<" not exist!!"<<endl;
      return 1;
    }

  for(;;)
    {
      if(fscanf(fin,"%lf %lf %lf",&lon,&lat,&temp)==EOF) break;
      if(lon<x0)x0=lon;
      if(lon>x1)x1=lon;
      if(lat<y0)y0=lat;
      if(lat>y1)y1=lat;
    }
  fclose(fin);
x0=floor(x0*10)/10.; 
y0=floor(y0*10)/10.;
x1=ceil(x1*10)/10.; 
y1=ceil(y1*10)/10.;
npts_x=int((x1-x0)/dx+1);
npts_y=int((y1-y0)/dy+1);

//printf("%f,%f,%f,%f,%d,%d",x0,x1,y0,y1,npts_x,npts_y);
//return 1;

  //fprintf(stderr,"Memory check!!\n");
  double tr_t[npts_x][npts_y],amp[npts_x][npts_y];
  double dx_km[npts_y],dy_km;
  //fprintf(stderr,"Memory enough!!\n");
  
//  x0=250;
//  y0=27;
//  x1=x0+(npts_x-1)*dx;
//  y1=y0+(npts_y-1)*dy;
  for(j=1;j<npts_y-1;j++)
    {
      lat_temp=y0+j*dy;
      lat_temp=atan(0.993277*tan(lat_temp/180*pi))*180/pi;
      dx_km[j]=radius*sin((90-lat_temp)/180*pi)*dx/180*pi;
    }
  dy_km=radius*dy/180*pi;
  sprintf(buff1,"%s.ph.txt_v1.HD",arg[1]);
  if((fin=fopen(buff1,"r"))==NULL)
    {
      cout<<buff1<<" not exist!!"<<endl;
      return 1;
    }
  sprintf(buff1,"%s_am.txt_v1.HD",arg[1]);
  if((fin_amp=fopen(buff1,"r"))==NULL)
    {
      cout<<buff1<<" not exist!!"<<endl;
      return 1;
    }

  for(i=0;i<npts_x;i++)
    {
      for(j=0;j<npts_y;j++)
	{
	 tr_t[i][j]=0;
	 amp[i][j]=0;
	}
    }
  
  for(;;)
    {
      if(fscanf(fin,"%lf %lf %lf",&lon,&lat,&temp)==EOF) break;
      if(lon>x1+0.01||lon<x0-0.01|lat>y1+0.01||lat<y0-0.01)
	continue;
      i=int((lon-x0)/dx+0.1);
      j=int((lat-y0)/dy+0.1);
      
      tr_t[i][j]=temp;
    }
  fclose(fin);
    
  for(;;)
    {
      if(fscanf(fin_amp,"%lf %lf %lf",&lon,&lat,&temp)==EOF) break;
      if(lon>x1+0.01||lon<x0-0.01|lat>y1+0.01||lat<y0-0.01)
        continue;
      i=int((lon-x0)/dx+0.1);
      j=int((lat-y0)/dy+0.1);

      amp[i][j]=temp;
    }
  fclose(fin_amp);


  //      double temp1,temp2;
  marker_nn=0;
  int marker;
  sprintf(buff1,"%s.ph.txt_v1",arg[1]);
  if((fin3=fopen(buff1,"r"))==NULL)
    {
      cout<<buff1<<" not exist!!"<<endl;
      return 1;
    } 
  //  sprintf(buff1,"%s.ph.txt_v2",arg[1]);
  //fout=fopen(buff1,"w");
  double arr[3000][5];
  int ii,jj,ist,ist_old, nam=0, nph=0, nvel=0;
  ist=0;
  ist_old=0;
  fgets(evbuff, 300, fin3);
  for(;;)
    {
      if(fscanf(fin3,"%lf %lf %lf %lf %lf %s",&lon2,&lat2,&temp2,&temp3,&temp4,stnm)==EOF) 
	{
	  fclose(fin3);
	  break;
	}
      ist_old++;
      if(temp2<2*period)
	{
	  arr[ist][0]=lon2;
	  arr[ist][1]=lat2;
	  arr[ist][2]=temp2;
	  arr[ist][3]=temp3;
	  arr[ist][4]=temp4;
	  //  fprintf(fout,"%lf %lf %lf %lf %lf\n",lon2,lat2,temp2,temp3,temp4);
	  ist++;
	  continue;
	}
      i=int((lon2-x0)/dx+0.5);
      j=int((lat2-y0)/dy+0.5);
      if(i<0||j<0||i>npts_x||j>npts_y)
	{
	  ist_old--;
	  //cout<<lon2<<" "<<lat2<<" "<<temp<<" out of range"<<endl;
	  continue;
	}
      temp=(tr_t[i+2][j]/-12.0+tr_t[i+1][j]*4.0/3.0+tr_t[i][j]*-5.0/2.0+tr_t[i-1][j]*4.0/3.0+tr_t[i-2][j]/-12.0)/dx_km[j]/dx_km[j];
      temp1=(tr_t[i][j+2]/-12.0+tr_t[i][j+1]*4.0/3.0+tr_t[i][j]*-5.0/2.0+tr_t[i][j-1]*4.0/3.0+tr_t[i][j-2]/-12.0)/dy_km/dy_km;
      temp=temp+temp1;
      if(temp>atof(arg[4])||temp<-(atof(arg[4])))
	{
	  //cout<<lon2<<" "<<lat2<<" curvature "<<temp<<endl;
          nph++;
	  continue;
	}
      temp=(tr_t[i+1][j]-tr_t[i-1][j])/2.0/dx_km[j];
      temp1=(tr_t[i][j+1]-tr_t[i][j-1])/2.0/dy_km;
      if(temp1==0)
	{
	  temp1=0.00001;
	}
      temp=sqrt(temp1*temp1+temp*temp);
      if(temp>1/2.0||temp<1/6.0)
	{
	  //cout<<lon2<<" "<<lat2<<" vel "<<temp<<endl;
          nvel++;
	  continue;
	}
      temp=(amp[i+2][j]/-12.0+amp[i+1][j]*4.0/3.0+amp[i][j]*-5.0/2.0+amp[i-1][j]*4.0/3.0+amp[i-2][j]/-12.0)/dx_km[j]/dx_km[j];
      temp1=(amp[i][j+2]/-12.0+amp[i][j+1]*4.0/3.0+amp[i][j]*-5.0/2.0+amp[i][j-1]*4.0/3.0+amp[i][j-2]/-12.0)/dy_km/dy_km;
      temp=fabs((temp+temp1)/amp[i][j]*period*period/4/pi/pi);
      if(temp>atof(arg[5]))
	{
	  //cout<<"amp_curvature "<<lon2<<" "<<lat2<<" "<<temp<<endl;
          nam++;
	  continue;
	}
      //marker=0;
      //for(ii=-3;ii<=3;ii+=3)
      //{
      //  for(jj=-3;jj<=3;jj+=3)
      //    {
      //      temp=(amp[i+ii+1][j+jj]-amp[i+ii-1][j+jj])/2.0/dx_km[j+jj];
      //      temp1=(amp[i+ii][j+jj+1]-amp[i+ii][j+jj-1])/2.0/dy_km;
      //      temp=sqrt(temp1*temp1+temp*temp)/amp[i+ii][j+jj]*period;
      //      if(temp>0.15)
      //	{
      //	  marker=1;
      //	  break;
      //	}
      //    }
      //  if(marker==1)
      //    break;
      //}
      //if(marker==1)
      //{
      //  cout<<lon2<<" "<<lat2<<" amp gradient > 0.15 "<<temp<<endl;
      //  continue;
      //}
      //      
      //cout<<"amp_curvature "<<lon2<<" "<<lat2<<" "<<temp<<endl;
      arr[ist][0]=lon2;
      arr[ist][1]=lat2;
      arr[ist][2]=temp2;
      arr[ist][3]=temp3;
      arr[ist][4]=temp4;
      //	  fprintf(fout,"%lf %lf %lf %lf %lf\n",lon2,lat2,temp2,temp3,temp4);
      ist++;
    }
  
  if(ist*1.0/ist_old<0.4)
    {
      cout<<"No result: "<<nam<<"(for amp) + "<<nph<<"(for ph) + "<<nvel<<"(for vel) out of "<<ist_old<<" stations removed."<<endl;
    }
  else
    if(ist<atof(arg[3]))
      {
	cout<<"No result: Only "<<ist<<" stations left!!"<<endl;
      }
    else
      {
	sprintf(buff1,"%s.ph.txt_v2",arg[1]);
	fout=fopen(buff1,"w");
        //fprintf(fout,"%s\n",evbuff);
	for(i=0;i<ist;i++)
	  {
	    fprintf(fout,"%lf %lf %lf %lf %lf\n",arr[i][0],arr[i][1],arr[i][2],arr[i][3],arr[i][4]); 
	  }
	fclose(fout);
        cout<<"Result OK."<<endl;
      }
  return 0;
}
