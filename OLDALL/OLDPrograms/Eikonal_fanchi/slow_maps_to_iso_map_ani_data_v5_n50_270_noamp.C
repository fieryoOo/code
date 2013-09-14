#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

int main(int na, char *arg[])
{
  if(na!=6)
    {
      cout<<"usage:travel_time_to_velocity_map station_morgan.lst min max N_bin out_name"<<endl;
      return 0;
    }
  FILE *ff,*fin,*fout,*file1,*file_iso,*file_ani;
  int i,j,k;
  int npts_x,npts_y;
  char buff1[300],sta1[10],name_iso[100],name_ani[100],name_ani_n[100];
  double lat,lon,t_lat,t_lon,radius,pi,sta1_lon,sta1_lat;
  int t_i,t_j,nsta;
  int ii,jj,kk,kkk,min_n;
  int marker_i,marker_j;
  double min,max,d_bin;
  int N_bin;
  min=atof(arg[2]);
  max=atof(arg[3]);
  N_bin=atoi(arg[4]);
  sprintf(name_iso,"%s.iso",arg[5]);
  sprintf(name_ani,"%s.ani",arg[5]);
  sprintf(name_ani_n,"%s_ani_n",arg[5]);  
  //  cout<<min<<" "<<max<<" "<<N_bin<<endl;
  double hist[N_bin];
  double slow_sum1[N_bin];
  double slow_un[N_bin];
  double azi_w2[N_bin], azi_std[N_bin];
  double azi_weight_sum[N_bin];
  //  for(i=0;i<N_bin;i++)
  //{
  //  hist[i]=0;
  //  slow_sum[i]=0;
  //  slow_un[i]=0;
  //}
  d_bin=(max-min)/N_bin;



  radius=6371.1391285;
  pi=4.0*atan(1.0);
  double dx,dy,x0,y0,x1,y1,temp,lat_temp,temp2,temp3,temp4,trash1,trash2;
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
  double slow[npts_x][npts_y][1000],slow_weight[npts_x][npts_y][1000];
  double azi[npts_x][npts_y][1000];
  double weight[1000];
  double weight_sum;
  int n[npts_x][npts_y];
  double slow_sum[npts_x][npts_y],slow_std[npts_x][npts_y];
  // double dx_km[npts_y],dy_km;
  
  fprintf(stderr,"Memory enough!!\n");
  
  dx=0.1;//degree
  dy=0.1;//degree
//  x0=250;
//  y0=27;
//  x1=x0+(npts_x-1)*dx;
//  y1=y0+(npts_y-1)*dy;
  //for(j=1;j<npts_y-1;j++)
  //  {
  //  lat_temp=y0+j*dy;
  //  lat_temp=atan(0.993277*tan(lat_temp/180*pi))*180/pi;
  //  dx_km[j]=radius*sin((90-lat_temp)/180*pi)*dx/180*pi;
  //}
  //dy_km=radius*dy/180*pi;
  for(i=0;i<npts_x;i++)
    {
      for(j=0;j<npts_y;j++)
	{
	  slow_sum[i][j]=0;
	  n[i][j]=0;
	}
    }
//  char event_name[300];
//  int marker_i,marker_j;
  marker_i=int((256.2-x0+0.001)/dx);
  marker_j=int((33.0-y0+0.001)/dy);
  file1=fopen(arg[1],"r");
  nsta=0;
  for(;;)
    {
      if(fscanf(file1,"%s",&event_name)==EOF)
        break;
      nsta++;
      sprintf(buff1,"slow_azi_%s.txt.HD",event_name);
      if((fin=fopen(buff1,"r"))==NULL)
       {
	 cout<<buff1<<" not exist!!"<<endl;
	 return 1;
       }
      for(;;)
       {
	 if(fscanf(fin,"%lf %lf %lf %lf",&lon,&lat,&temp,&temp2)==EOF) break;
	 if(lon>x1+0.01||lon<x0-0.01|lat>y1+0.01||lat<y0-0.01)
	   continue;
	 i=int((lon-x0)/dx+0.1);
	 j=int((lat-y0)/dy+0.1);
	 if(temp<0.5&&temp>0.1)
	   {
	     slow[i][j][n[i][j]]=temp;
	     //slow[i][j][n[i][j]]=sqrt(temp*temp-temp4);
	     azi[i][j][n[i][j]]=temp2;
	     //     slow_weight[i][j][n[i][j]]=fabs(0.5*temp4/slow[i][j][n[i][j]])+0.01;
	     // slow_weight[i][j][n[i][j]]=1/slow_weight[i][j][n[i][j]]/slow_weight[i][j][n[i][j]];
	     
	     //	     if(i==marker_i&&j==marker_j)
	     //fprintf(stderr,"%lf %lf %lf\n",1/temp,slow_weight[i][j][n[i][j]],temp4);
	     slow_weight[i][j][n[i][j]]=1;
	     n[i][j]++;
	     //	     slow_sum[i][j]+=temp;
	   }
       }
      fclose(fin);
    }
  fclose(file1);
  file_iso=fopen(name_iso,"w");
  nsta=20;
  double w2;
  double temp_slow_sum;
  int temp_n;
  for(i=0;i<npts_x;i++)
    {
      for(j=0;j<npts_y;j++)
	{
	  if(n[i][j]<0.5*nsta)
	    {
	      fprintf(file_iso,"%lf %lf 0 9999 %d\n",x0+i*dx,y0+j*dy,n[i][j]);
	      continue;
	    }
	  w2=0;
	  weight_sum=0;
	  for(k=0;k<n[i][j];k++)
	    {
	      weight[k]=0;
	      for(kk=0;kk<n[i][j];kk++)
		{
		  if(fabs(azi[i][j][kk]-azi[i][j][k])<25)
		    weight[k]++;
		}
	      weight[k]=1/weight[k]*slow_weight[i][j][k];
	      weight_sum+=weight[k];
	      //if(i==marker_i&&j==marker_j)
	      //fprintf(stderr,"%lf %lf %lf\n",1/slow[i][j][k],weight[k],slow_weight[i][j][k]);
	    }
	  for(k=0;k<n[i][j];k++)
	    {
	      weight[k]=weight[k]/weight_sum;
	      slow_sum[i][j]+=weight[k]*slow[i][j][k];
	      w2+=weight[k]*weight[k];
	    }
	  //	  if(n[i][j]>=0.5*nsta)
	  //{
	  //  
	  //  slow_sum[i][j]=slow_sum[i][j]/n[i][j];
	  temp_slow_sum=slow_sum[i][j];
	  temp=0;
	  for(k=0;k<n[i][j];k++)
	    {
	      temp+=weight[k]*(slow[i][j][k]-slow_sum[i][j])*(slow[i][j][k]-slow_sum[i][j]);
	    }
	  slow_std[i][j]=sqrt(temp/(1-w2));
	  // if(i==marker_i&&j==marker_j)
	  //fprintf(stderr,"%lf\n",slow_std[i][j]);
	  w2=0;
	  weight_sum=0;
	  slow_sum[i][j]=0;
	  temp_n=0;
	  for(k=0;k<n[i][j];k++)
	    {
	      if(fabs(slow[i][j][k]-temp_slow_sum)>2.0*slow_std[i][j])
		continue;
	      weight_sum+=weight[k];
	      temp_n++;
	    }
	  //	   if(i==marker_i&&j==marker_j)
	  // fprintf(stderr,"%d %d\n",n[i][j],temp_n);
	  for(k=0;k<n[i][j];k++)
	    {
	      if(fabs(slow[i][j][k]-temp_slow_sum)>2.0*slow_std[i][j])
		continue;
	      weight[k]=weight[k]/weight_sum;
	      slow_sum[i][j]+=weight[k]*slow[i][j][k];
	      w2+=weight[k]*weight[k];
	    }
	  temp=0;
	  for(k=0;k<n[i][j];k++)
	    {
	      if(fabs(slow[i][j][k]-temp_slow_sum)>2.0*slow_std[i][j])
		continue;
	      temp+=weight[k]*(slow[i][j][k]-slow_sum[i][j])*(slow[i][j][k]-slow_sum[i][j]);
	    }
	  slow_std[i][j]=sqrt(temp/(1-w2));

	  temp=slow_std[i][j]*sqrt(w2)/slow_sum[i][j]/slow_sum[i][j];
	  //	  temp=sqrt(temp/(n[i][j]-1)/n[i][j])/slow_sum[i][j]/slow_sum[i][j];
	  // cout<<x0+i*dx<<" "<<y0+j*dy<<" "<<1/slow_sum[i][j]<<" "<<temp<<" "<<n[i][j]<<endl;
	  fprintf(file_iso,"%lf %lf %lf %lf %d\n",x0+i*dx,y0+j*dy,1/slow_sum[i][j],temp,temp_n);
	  //}
	  //else
	  //{
	  //  //    cout<<x0+i*dx<<" "<<y0+j*dy<<" 0 "<<" 9999 "<<n[i][j]<<endl;
	  //  fprintf(file_iso,"%lf %lf 0 9999 %d\n",x0+i*dx,y0+j*dy,n[i][j]);
	  //}
	}
    }
  fclose(file_iso);
  file_ani=fopen(name_ani,"w");
  fout=fopen(name_ani_n,"w");  


  double test_lon,test_lat;
  //test_lon=240.6;
  //test_lat=36.6;
  
  for(i=0;i<npts_x;i++)
    {
      for(j=0;j<npts_y;j++)
        {
	  //	  fprintf(stderr,"gill!! %d %d\n",i,j);
	  //  i=int((test_lon-x0)/dx+0.1);
	  //j=int((test_lat-y0)/dy+0.1);
	  //cout<<i<<" "<<j<<endl;
	  //fprintf(file_ani,"\n>\n%lf %lf\n",x0+i*dx,y0+j*dy);
	  for(k=0;k<N_bin;k++)
	    {
	      hist[k]=0;
	      slow_sum1[k]=0;
	      slow_un[k]=0;
	      azi_weight_sum[k]=0;
	      azi_w2[k]=0;
	      azi_std[k]=0;
	    }
	  if(i-3<0||i+3>=npts_x||j-3<0||j+3>=npts_y)
	    continue;
	  kkk=0;
	  //	  min_n=999999999;
	  for(ii=i-3;ii<=i+3;ii+=3)
	    {
	      for(jj=j-3;jj<=j+3;jj+=3)
		{
		  if(n[ii][jj]<nsta*0.5)
		    continue;
		  kkk+=n[ii][jj];
		  //  if(n[ii][jj]<min_n)
		  //min_n=n[ii][jj];
		}
	    }
	  fprintf(fout,"%lf %lf %d\n",x0+i*dx,y0+j*dy,kkk);
	  //fprintf(stderr,"%lf %lf %d\n",x0+i*dx,y0+j*dy,kkk);
	  if(kkk<9*nsta*0.5||n[i][j]<nsta*0.5)
	    continue;
	  for(ii=i-3;ii<=i+3;ii+=3)
	    {
	      for(jj=j-3;jj<=j+3;jj+=3)
		{
		  if(n[ii][jj]<nsta*0.5)
		    continue;
		  for(k=0;k<n[ii][jj];k++)
		    {
		      if(azi[ii][jj][k]>max||azi[ii][jj][k]<min)
			{
			  fprintf(stderr,"out of range!!");
			  return 1;
			}
		      //		      fprintf(stderr,"oliver!! %d %d\n",ii,jj);
		      hist[int((azi[ii][jj][k]-min)/d_bin)]++;
		      azi_weight_sum[int((azi[ii][jj][k]-min)/d_bin)]+=slow_weight[ii][jj][k];
		      azi_w2[int((azi[ii][jj][k]-min)/d_bin)]+=slow_weight[ii][jj][k]*slow_weight[ii][jj][k];
		      slow_sum1[int((azi[ii][jj][k]-min)/d_bin)]+=(slow[ii][jj][k]-slow_sum[ii][jj])*slow_weight[ii][jj][k];
		      //slow_un[int((azi[ii][jj][k]-min)/d_bin)]+=(slow[ii][jj][k]-slow_sum[ii][jj])*(slow[ii][jj][k]-slow_sum[ii][jj]);
		    }
		}
	    }
	  //	  fprintf(stderr,"tobal!! %d %d\n",i,j);
	  kk=0;
	  
	  for(k=0;k<N_bin;k++)
	    {
	      if(hist[k]>=10)
		{
		  kk++;
		  azi_w2[k]=azi_w2[k]/azi_weight_sum[k]/azi_weight_sum[k];
		  slow_sum1[k]=slow_sum1[k]/azi_weight_sum[k];
		}
	    }
	  //	  fprintf(stderr,"tobal_1!! %d %d\n",i,j);

	  for(ii=i-3;ii<=i+3;ii+=3)
            {
              for(jj=j-3;jj<=j+3;jj+=3)
                {
                  if(n[ii][jj]<nsta*0.5)
                    continue;
                  for(k=0;k<n[ii][jj];k++)
                    {
                      if(azi[ii][jj][k]>max||azi[ii][jj][k]<min)
                        {
                          fprintf(stderr,"out of range!!");
                          return 1;
                        }
		      azi_std[int((azi[ii][jj][k]-min)/d_bin)]+=slow_weight[ii][jj][k]/azi_weight_sum[int((azi[ii][jj][k]-min)/d_bin)]*(slow[ii][jj][k]-slow_sum[ii][jj]-slow_sum1[int((azi[ii][jj][k]-min)/d_bin)])*(slow[ii][jj][k]-slow_sum[ii][jj]-slow_sum1[int((azi[ii][jj][k]-min)/d_bin)]);
                    }
                }
            }
	  //	  fprintf(stderr,"tobal_2!! %d %d\n",i,j);

	  
	  fprintf(file_ani,"%lf %lf %d\n",x0+i*dx,y0+j*dy,kk);
	  for(k=0;k<N_bin;k++)
	    {
	      if(hist[k]>=10)
		{
		  //		  slow_sum1[k]=slow_sum1[k]/hist[k];
		  //slow_sum1[k]=slow_sum1[k]/azi_weight_sum[k];
		  //		  slow_un[k]=slow_un[k]-slow_sum1[k]*slow_sum1[k]*hist[k];
		  //slow_un[k]=sqrt(slow_un[k]/(hist[k]-1)/hist[k])/(slow_sum[i][j]+slow_sum1[k])/(slow_sum[i][j]+slow_sum1[k]); //uncertainty of vel not slow
		  //slow_un[k]=slow_std[i][j]/sqrt(double(hist[k]));
		  //slow_un[k]=slow_un[k]/(slow_sum[i][j]+slow_sum1[k])/(slow_sum[i][j]+slow_sum1[k]);//uncertainty of vel not slow
		  //slow_un[k]=slow_std[i][j];
		  
		  azi_std[k]=sqrt(azi_std[k]/(1-azi_w2[k]));
		  temp=azi_std[k]*sqrt(azi_w2[k])/slow_sum[i][j]/slow_sum[i][j];
		  //slow_un[k]=9999999;
		  //  to here             ---------
		  //cout<<min+(0.5+k)*d_bin<<" "<<1/slow_sum1[k]<<" "<<slow_un[k]<<endl;
		  //fprintf(file_ani,"%lf %lf %lf\n",min+(0.5+k)*d_bin,1/(slow_sum[i][j]+slow_sum1[k]),slow_un[k]);
		  fprintf(file_ani,"%lf %lf %lf\n",min+(0.5+k)*d_bin,1/(slow_sum[i][j]+slow_sum1[k]),temp);
		}
	    }
	  //	  fprintf(stderr,"done!!\n");
	  //return 0;
	}
    }
  fclose(file_ani);
  fclose(fout);
  return 0;
}
