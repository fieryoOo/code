#define MAIN
#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>


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
      std::cout<<"warning cos(theta)>1 and correct to 1!!"<<temp<<std::endl;
      temp=1;
    }
  if(temp<-1)
    {
      std::cout<<"warning cos(theta)<-1 and correct to -1!!"<<temp<<std::endl;
      temp=-1;
    }
  theta=fabs(acos(temp));
  return theta*radius;
}


int main(int na, char *arg[])
{ 
  if(na!=4)
    {
      std::cout<<"usage:correct_2pi event_file_name period least_sta_#"<<std::endl;
      return 0;
    }
  FILE *ff,*f2;
  double dis,dis_min;
  double per;
  per=atof(arg[2]);
  std::ifstream f1(arg[1]);
  if( ! f1 ) {
     std::cerr<<"Cannot read from file "<<arg[1]<<std::endl;
     exit(1);
  }
  int N=3000;
  double lon[N],lat[N],time[N],vel[N],amp[N];
  int YesNo[N];
  int i,j,k,ii,jj,iev,ist;
  // cout<<"test"<<endl;
  i=0;
  for( std::string line; std::getline(f1, line); ) {
     if( sscanf(line.c_str(),"%lf %lf %lf %lf %lf",&lon[i],&lat[i],&time[i],&vel[i],&amp[i]) != 5 ) continue;
     YesNo[i]=0; i++;
  }
  ist=i;

  //cout<<ist<<endl;
  f1.close();
  if(ist<atof(arg[3]))
    {
      printf("Only %d stations\n",ist);
      exit(1);
    }
  /* sort by distance to the center */
  int order[ist];
  double fact[ist];
  //double slon=245.133, slat=41.153;
  double slon=260, slat=30;
  for(i=0;i<ist;i++)
    {
      //      fact[i]=lat[i]-lon[i];
      fact[i]=(lat[i]-slat)*(lat[i]-slat)+(lon[i]-slon)*(lon[i]-slon);
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

  int mark,temp;
  char buff[300];
  sprintf(buff,"%s_v1",arg[1]);
  //ff=fopen(buff,"w");
  // fprintf(ff,"%lf %lf %lf %lf %lf\n",lon[order[0]],lat[order[0]],time[order[0]],vel[order[0]],amp[order[0]]);
  int ino=0;
  //  for(i=0;i<ist;i++)
  //cout<<i<<" "<<lon[order[i]]<<" "<<lat[order[i]]<<endl;
  
//for(int i=0; i<ist; i++) if(vel[i]<2.3)std::cerr<<"vel<2.3: "<<vel[i]<<std::endl;
  for(i=1;i<ist;i++)
    {
//std::cerr<<"233: "<<lon[233]<<" "<<lat[233]<<" "<<vel[233]<<" "<<amp[233]<<std::endl;
      /* search for a closerst station among those corrected */
      dis_min=999999999;
      for(j=0;j<i;j++)
	{
	  //if(j==137)
	    //	    fprintf(stderr,"%d %d %lf %lf\n",i,j,lon[order[j]],lat[order[j]]);
	  //	  cout<<lat[order[i]]<<" "<<lon[order[i]]<<" "<<lat[order[j]]<<" "<<lon[order[j]]<<endl;
	  dis=get_dist(lat[order[i]],lon[order[i]],lat[order[j]],lon[order[j]]);
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
      //      if(i==326)
      //fprintf(stderr,"%d %d %lf %lf %lf %lf\n",i,mark);
      //      cout<<">"<<endl;
      //cout<<lon[order[i]]<<" "<<lat[order[i]]<<endl;
     //cout<<lon[order[mark]]<<" "<<lat[order[mark]]<<endl;
      if(fabs(amp[order[i]]-amp[order[mark]])>10*amp[order[i]]||fabs(amp[order[i]]-amp[order[mark]])>10*amp[order[mark]])
	{
/*
std::cerr<<"order[i]="<<order[i]<<"  order[mark]="<<order[mark]<<std::endl;
std::cerr<<" ("<<lon[order[i]]<<", "<<lat[order[i]]<<") - ("<<lon[order[mark]]<<", "<<lat[order[mark]]<<") "<<std::endl;
std::cerr<<" amp1="<<amp[order[i]]<<"  amp2="<<amp[order[mark]]<<std::endl;
std::cerr<<" vel1="<<vel[order[i]]<<"  vel2="<<vel[order[mark]]<<std::endl;
std::cerr<<"before: amp="<<amp[order[i]]<<" vel="<<vel[order[i]]<<std::endl;
*/
	  //	  if(i==137)
	  // fprintf(stderr,"%d %d\n",i,mark);
	  std::cout<<dis<<" "<<amp[order[i]]<<" "<<amp[order[mark]]<<std::endl;
	  YesNo[order[i]]=999;
	  ino++;
	  lon[order[i]]=lon[order[mark]];
	  lat[order[i]]=lat[order[mark]];
	  time[order[i]]=time[order[mark]];
	  vel[order[i]]=vel[order[mark]];
	  amp[order[i]]=amp[order[mark]];
//std::cerr<<"after: amp="<<amp[order[i]]<<" vel="<<vel[order[i]]<<std::endl;
	  continue;
	}
      dis=vel[order[i]]*time[order[i]];
      temp=dis/vel[order[mark]];
      //std::cerr<<"per="<<per<<" dis="<<dis<<" vel="<<vel[order[mark]]<<" Tpred="<<temp<<" Tobs="<<time[order[i]]<<std::endl;
      for(;;)
	{
	  if(time[order[i]]-temp>per/2)
	    time[order[i]]-=per;
	  else
	    if(temp-time[order[i]]>per/2)
	      time[order[i]]+=per;
	    else
	      break;
	}
      if(fabs(time[order[i]]-temp)>6)
	{
	  //	  if(i==137)
	  //fprintf(stderr,"%d %d %lf %lf\n",i,mark,lon[order[mark]],lat[order[mark]]);
	  std::cout<<"misft > 6 sec"<<std::endl;
	  YesNo[order[i]]=999;
	  ino++;
	  lon[order[i]]=lon[order[mark]];
	  lat[order[i]]=lat[order[mark]];
	  time[order[i]]=time[order[mark]];
	  vel[order[i]]=vel[order[mark]];
	  amp[order[i]]=amp[order[mark]];
	  continue;
	}
std::cerr<<"vel0="<<vel[order[i]];
      vel[order[i]]=dis/time[order[i]];
std::cerr<<" velc="<<vel[order[i]]<<std::endl;
      //fprintf(ff,"%lf %lf %lf %lf %lf\n",lon[order[i]],lat[order[i]],time[order[i]],vel[order[i]],amp[order[i]]);
      //      fprintf(ff,"%lf %lf %lf %lf %lf %lf %lf %lf\n",lon[order[i]],lat[order[i]],time[order[i]],vel[order[i]],amp[order[i]],lon[order[mark]],lat[order[mark]],vel[order[mark]]);
      //cout<<lon[order[i]]<<" "<<lat[order[i]]<<" "<<fact[order[i]]<<endl;
    }
//for(int i=0; i<ist; i++) if(vel[i]<2.3)std::cerr<<"vel<2.3: "<<vel[i]<<std::endl;
  //  cout<<ino<<" "<<ist<<" "<<ino/ist-0.9<<endl;
  if(ino*1.0/ist>0.5||ist-ino<atof(arg[3]))
    {
      printf("too many stations removed!! %d %d\n",ist,ino);
      exit(1);
    }
  ff=fopen(buff,"w");
  
  for(i=0;i<ist;i++)
    {
      if(YesNo[order[i]]==0)
	fprintf(ff,"%lf %lf %lf %lf %lf\n",lon[order[i]],lat[order[i]],time[order[i]],vel[order[i]],amp[order[i]]);	
    }
  
  fclose(ff);
  return 0;
}
