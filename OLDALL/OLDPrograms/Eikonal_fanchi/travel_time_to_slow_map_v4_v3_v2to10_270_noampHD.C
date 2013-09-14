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
      //cout<<lat1<<" "<<lon1<<" "<<lat2<<" "<<lon2<<endl;
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
  if(na!=3)
    {
      cout<<"usage:travel_time_to_slow_map event.lst period"<<endl;
      return 0;
    }
  FILE *ff,*fin,*fin2,*fin3,*fout,*file1,*fin_amp;
  int i,j,ii,jj,iii,jjj;
  int npts_x,npts_y;
  char buff1[300],sta1[10];
  double lat,lon,lat2,lon2,t_lat,t_lon,radius,pi;//sta1_lon,sta1_lat;
  double sta1_lon[1000],sta1_lat[1000],t0[1000];
  int t_i,t_j;
  int marker_nn,marker_EN[2][2],marker_E,marker_N;
  double period,dist;
  int nn,nsta;
  period=atof(arg[2]);
  
  radius=6371.1391285;
  pi=4.0*atan(1.0);
  double w2;
  w2=(2*pi/period)*(2*pi/period);
  double dx,dy,x0,y0,x1,y1,temp,temp1,temp2,lat_temp, temp3, temp4,temp_amp;
//  npts_x=63;
//  npts_y=61;
  dx=0.1;//degree
  dy=0.1;//degree
 
 x0=999;
 x1=-999;
 y0=999;
 y1=-999;

  char event_name[300];

  file1=fopen(arg[1],"r");
  for(;;)
  {
      if(fscanf(file1,"%s",&event_name)==EOF)
        break;

  sprintf(buff1,"%s.ph.txt_v1.HD",event_name);
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
  }
fclose(file1);
x0=floor(x0*10)/10.;
y0=floor(y0*10)/10.;
x1=ceil(x1*10)/10.;
y1=ceil(y1*10)/10.;
npts_x=int((x1-x0)/dx+1);
npts_y=int((y1-y0)/dy+1);

//printf("%f,%f,%f,%f,%d,%d",x0,x1,y0,y1,npts_x,npts_y);
//return 1;

  fprintf(stderr,"Memory check!!\n");
  double tr_t[npts_x][npts_y],amp[npts_x][npts_y],amp_grad[npts_x][npts_y];
  double dx_km[npts_y],dy_km;
  double slow_all[npts_x][npts_y][4];

  fprintf(stderr,"Memory enough!!\n");
  
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
  
//  char event_name[300];
  
  file1=fopen(arg[1],"r");
  for(;;)
    {
      if(fscanf(file1,"%s",&event_name)==EOF)
	break;
      sprintf(buff1,"slow_azi_%s.txt.HD",event_name);
      if(access(buff1, F_OK) == 0)
	{
	  fprintf(stderr,"%s exist!! skip\n",buff1);
	  continue;
	}
      sprintf(buff1,"%s.ph.txt_v2.HD",event_name);
      if((fin=fopen(buff1,"r"))==NULL)
       {
	 cout<<buff1<<" not exist!!"<<endl;
	 return 1;
       }
      sprintf(buff1,"%s.ph.txt_v2.HD_0.2",event_name);
      if((fin2=fopen(buff1,"r"))==NULL)
	{
	  cout<<buff1<<" not exist!!"<<endl;
	  return 1;
	}
      //      sprintf(buff1,"%s_am.txt_v2.HD",event_name);
      //if((fin_amp=fopen(buff1,"r"))==NULL)
      //{
      //  cout<<buff1<<" not exist!!"<<endl;
      //  return 1;
      //}

      sprintf(buff1,"slow_azi_%s.txt.HD",event_name);
      fout=fopen(buff1,"w");
      for(i=0;i<npts_x;i++)
	{
	  for(j=0;j<npts_y;j++)
	    {
	      tr_t[i][j]=0;
	      amp[i][j]=0;
	      amp_grad[i][j]=9999;
	    }
	}
      sprintf(buff1,"%s.ph.txt_v2",event_name);
      if((fin3=fopen(buff1,"r"))==NULL)
	{
	  cout<<buff1<<" not exist!!"<<endl;
	  return 1;
	}    
      nsta=0;
      for(;;)
	{
	  if(fscanf(fin3,"%lf %lf %lf %lf %lf",&sta1_lon[nsta],&sta1_lat[nsta],&temp2,&temp3,&temp4)==EOF) 
	    break;
	  nsta++;
	  if(nsta>=1000)
	    {
	      fprintf(stderr,"too many stations!!\n");
	      return 0;  
	    }
	}
      fclose(fin3);
      for(;;)
	{
	  if(fscanf(fin,"%lf %lf %lf",&lon,&lat,&temp)==EOF) break;
	  if(fscanf(fin2,"%lf %lf %lf",&lon2,&lat2,&temp2)==EOF) break;
	  if(lon!=lon2||lat!=lat2)
	    {
	      fprintf(stderr,"HD and HD_0.2 files not compatiable!!\n");
	      return 0;
	    }
	  if(lon>x1+0.01||lon<x0-0.01|lat>y1+0.01||lat<y0-0.01)
	    {
	      //printf("grid point pass the boundary!!\n"); 
	      continue;
	    }
	  i=int((lon-x0)/dx+0.1);
	  j=int((lat-y0)/dy+0.1);
	  if(temp<temp2-1||temp>temp2+1||temp<2*period)
	    {
	      tr_t[i][j]=0;	     
	      continue;
	    }
	  marker_nn=3;
	  marker_EN[0][0]=0;
	  marker_EN[0][1]=0;
	  marker_EN[1][0]=0;
          marker_EN[1][1]=0;
	  // fin3=fopen(buff1,"r");
	  for(ii=0;ii<=nsta;ii++)
	    {
	      //	      fprintf(stderr,"test!!\n");
	      //	      if(fscanf(fin3,"%lf %lf %lf %lf %lf",&lon2,&lat2,&temp2,&temp3,&temp4)==EOF) 
	      //{
	      if(ii==nsta)
		{
		  temp=0;
		  break;
		}
	      //  fclose(fin3);
	      //  break;
	      //}
	      //	      cout<<"test!! "<<lon2<<" "<<lat2<<endl;
	      if(sta1_lon[ii]-lon<0)
		marker_E=0;
	      else
		marker_E=1;
	      if(sta1_lat[ii]-lat<0)
                marker_N=0;
              else
                marker_N=1;
	      if(marker_EN[marker_E][marker_N]!=0)
		continue;
	      //cout<<lat<<" "<<lon<<" "<<lat2<<" "<<lon2<<endl;
	      dist=get_dist(lat,lon,sta1_lat[ii],sta1_lon[ii]);
	      if(dist<150)
		{
		  marker_nn--;
		  if(marker_nn==0)
		    {
		      //      fclose(fin3);
		      break;
		    }
		  marker_EN[marker_E][marker_N]++;
		}
	    }
	  tr_t[i][j]=temp;
	  //cout<<temp<<endl;
	}
      fclose(fin);
      fclose(fin2);
      //      fprintf(stderr,"test!!\n");
      //for(;;)
      //{
      //  if(fscanf(fin_amp,"%lf %lf %lf",&lon,&lat,&temp)==EOF) break;
      //  if(lon>x1+0.01||lon<x0-0.01|lat>y1+0.01||lat<y0-0.01)
      //    continue;
      //  i=int((lon-x0)/dx+0.1);
      //  j=int((lat-y0)/dy+0.1);
      //  
      //  amp[i][j]=temp;
      //}
      //fclose(fin_amp);
      //      for(i=1;i<npts_x-1;i++)
      //{
      //  for(j=1;j<npts_y-1;j++)
      //    {
      //      temp1=(amp[i+1][j]-amp[i-1][j])/2.0/dx_km[j];
      //      temp2=(amp[i][j+1]-amp[i][j-1])/2.0/dy_km;
      //      amp_grad[i][j]=sqrt(temp1*temp1+temp2*temp2)/amp[i][j]*period;
      //    }
      //}
      //      double temp1,temp2;
      for(i=4;i<npts_x-4;i++)
	{
	  for(j=4;j<npts_y-4;j++)
	    {
	      temp1=(tr_t[i+1][j]-tr_t[i-1][j])/2.0/dx_km[j];
	      temp2=(tr_t[i][j+1]-tr_t[i][j-1])/2.0/dy_km;
	      if(temp2==0)
		{
		  temp2=0.00001;
		}
	      temp=sqrt(temp1*temp1+temp2*temp2);
printf("for temp1: %lf-%lf/(2*%lf)\n",tr_t[i+1][j],tr_t[i-1][j],dx_km[j]);
printf("for temp2: %lf-%lf/(2*%lf)\n",tr_t[i][j+1],tr_t[i][j-1],dy_km);
printf("0.1<%lf<0.5?\n",temp);

	      if(temp>0.5||temp<0.1||tr_t[i+1][j]==0||tr_t[i-1][j]==0||tr_t[i][j+1]==0||tr_t[i][j-1]==0)
		{
		  slow_all[i][j][0]=0;
		  slow_all[i][j][1]=999;
		  slow_all[i][j][2]=999;
		  slow_all[i][j][3]=999;
		  //fprintf(fout,"%lf %lf 0 999 999 999\n",x0+i*dx,y0+j*dy);
		}
	      
	      else
		{
		  //		  temp3=(amp[i+2][j]-amp[i-2][j])/2.0/dx_km[j]/2;
		  //temp4=(amp[i][j+2]-amp[i][j-2])/2.0/dy_km/2;
		  //temp_amp=sqrt(temp1*temp1+temp2*temp2)/amp[i][j]*period;
		  //temp_amp=0;
		  //	  for(ii=-3;ii<=3;ii++)
		  //for(jj=-3;jj<=3;jj++)
		  //  {
		  //iii=i+ii;
		  //jjj=j+jj;
		  //if(iii<1||iii>npts_x-1||jjj<1||jjj>npts_y-1)
		  //  continue;
		  //if(amp_grad[iii][jjj]>temp_amp)
		  //  temp_amp=amp_grad[iii][jjj];
		  //		  temp3=(tr_t[i+2][j]/-12.0+tr_t[i+1][j]*4.0/3.0+tr_t[i][j]*-5.0/2.0+tr_t[i-1][j]*4.0/3.0+tr_t[i-2][j]/-12.0)/dx_km[j]/dx_km[j];
		  //temp4=(tr_t[i][j+2]/-12.0+tr_t[i][j+1]*4.0/3.0+tr_t[i][j]*-5.0/2.0+tr_t[i][j-1]*4.0/3.0+tr_t[i][j-2]/-12.0)/dy_km/dy_km;
		  //temp3=(amp[i+4][j]/-12.0+amp[i+2][j]*4.0/3.0+amp[i][j]*-5.0/2.0+amp[i-2][j]*4.0/3.0+amp[i-4][j]/-12.0)/dx_km[j]/dx_km[j]/4;
                  //temp4=(amp[i][j+4]/-12.0+amp[i][j+2]*4.0/3.0+amp[i][j]*-5.0/2.0+amp[i][j-2]*4.0/3.0+amp[i][j-4]/-12.0)/dy_km/dy_km/4;
		  //temp4=(temp3+temp4)/amp[i][j]/w2;
		  
		  if(temp1>0&&temp2>0)
		    {
		      slow_all[i][j][0]=temp;
		      slow_all[i][j][1]=atan(temp2/temp1)/pi*180;
		      //  slow_all[i][j][2]=temp_amp;
		      //slow_all[i][j][3]=temp4;
		      //		      fprintf(fout,"%lf %lf %lf %lf %lf %lf\n",x0+i*dx,y0+j*dy,temp,atan(temp2/temp1)/pi*180,temp_amp,temp4);
		    }
		  if(temp1>0&&temp2<=0)
		    {
		      slow_all[i][j][0]=temp;
		      slow_all[i][j][1]=360+atan(temp2/temp1)/pi*180;
		      //slow_all[i][j][2]=temp_amp;
		      //slow_all[i][j][3]=temp4;
		      //fprintf(fout,"%lf %lf %lf %lf %lf %lf\n",x0+i*dx,y0+j*dy,temp,360+atan(temp2/temp1)/pi*180,temp_amp,temp4);
		    }
		  if(temp1<=0)
		    { 
		      slow_all[i][j][0]=temp;
		      slow_all[i][j][1]=180+atan(temp2/temp1)/pi*180;
		      //slow_all[i][j][2]=temp_amp;
		      //slow_all[i][j][3]=temp4;
		      //		    fprintf(fout,"%lf %lf %lf %lf %lf %lf\n",x0+i*dx,y0+j*dy,temp,180+atan(temp2/temp1)/pi*180,temp_amp,temp4);
		    }		
		}
	    }
	}
      for(i=4;i<npts_x-4;i++)
	{
	  for(j=4;j<npts_y-4;j++)
	    {
	      if(slow_all[i][j][0]==0)
		{
		  fprintf(fout,"%lf %lf 0 999\n",x0+i*dx,y0+j*dy);
		  continue;
		}
	      nn=0;
	      for(ii=i-3;ii<=i+3;ii+=3)
		{
		  for(jj=j-3;jj<=j+3;jj+=3)
		    {
		      if(slow_all[ii][jj][0]!=0)
			nn++;
		    }
		}
	      if(nn>3)
		fprintf(fout,"%lf %lf %lf %lf\n",x0+i*dx,y0+j*dy,slow_all[i][j][0],slow_all[i][j][1]);
	      else
		fprintf(fout,"%lf %lf 0 999\n",x0+i*dx,y0+j*dy);
	    }
	}
      fclose(fout);
    }
  fclose(file1);
  return 0;
}
