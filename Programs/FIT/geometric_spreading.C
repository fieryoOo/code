#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;

int main (int argc, char *argv[])
{
  FILE *f;
  int i,j,k,ii,idata,nstep,step;
  double slp,slpmx,slpmn,slpstp,R,Rmx,Rmn,Rstp,itcpt,itcptmx,itcptmn,itcptstp,slpfit,itcptfit,Rfit,slpp[1000],Rr[1000],itcptt[1000],datax[5000],datay[5000],rsd[1000],dataxx,datayy,rsdmn,rsdfit,tempc;

  if(argc != 2)
    {
      printf("Usage: geometric_spreading.C [input file]\n");
      exit(-1);
    }

  if((f = fopen(argv[1],"r"))==NULL){
     fprintf(stderr,"cannot open data file %s\n",argv[1]);
     return 0;
    }
  for (i=0;;i++) {
     if (fscanf(f,"%lf %lf", &datax[i],&datay[i]) != 2) break;
    }
  idata=i;
  fclose(f);

  for(i=0;i<idata;i++)
     for(j=i;j<idata;j++)
        if(datax[i]>datax[j]){
           dataxx=datax[i];
           datayy=datay[i];
           datax[i]=datax[j];
           datay[i]=datay[j];
           datax[j]=dataxx;
           datay[j]=datayy;
          }

  slpmx=0.;slpmn=-0.01;
  Rmx=20000.;Rmn=0.;
  itcptmx=-999999999.;itcptmn=999999999.;
  for(i=0;i<idata;i++){
     itcpt=datay[i]-datax[i]*slpmx;
     if(itcptmn>itcpt)itcptmn=itcpt;
     itcpt=datay[i]-datax[i]*slpmn;
     if(itcptmx<itcpt)itcptmx=itcpt;
    }

  for(step=0;;step++){
//printf("slpmin: %f  slpmax: %f\n",slpmn,slpmx);
    slpstp=(slpmx-slpmn)/9;
    itcptstp=(itcptmx-itcptmn)/9;
    Rstp=(Rmx-Rmn)/9;
//printf("slpstp: %f  itcptstp: %f  Rstp: %f\n",slpstp,itcptstp,Rstp);
  
    rsdmn=999999999.;
    for(i=0;i<10;i++){
       slp=slpmn+i*slpstp;
       for(j=0;j<10;j++){
          itcpt=itcptmn+j*itcptstp;
          for(k=0;k<10;k++){
             R=Rmn+k*Rstp;
             rsd[i*100+j*10+k]=0;
             slpp[i*100+j*10+k]=slp;
             itcptt[i*100+j*10+k]=itcpt;
             Rr[i*100+j*10+k]=R;
             for(ii=0;ii<idata;ii++){
                rsd[i*100+j*10+k]+=(slp*datax[ii]+itcpt+0.5*log(datax[ii]/(datax[ii]+R))-datay[ii])*(slp*datax[ii]+itcpt+0.5*log(datax[ii]/(datax[ii]+R))-datay[ii]);
               }
             rsd[i*100+j*10+k]=sqrt(rsd[i*100+j*10+k]/(idata-1));
            }
         }
      }  
    for(i=0;i<5;i++)
       for(j=i;j<1000;j++)
          if(rsd[i]>rsd[j]){
             rsdfit=rsd[i];
             slpfit=slpp[i];
             itcptfit=itcptt[i];
             Rfit=Rr[i];
             rsd[i]=rsd[j];
             slpp[i]=slpp[j];
             itcptt[i]=itcptt[j];
             Rr[i]=Rr[j];
             rsd[j]=rsdfit;
             slpp[j]=slpfit;
             itcptt[j]=itcptfit;
             Rr[j]=Rfit;
            }
  
    tempc=slpmn;slpmn=slpmx;slpmx=tempc;
    tempc=itcptmn;itcptmn=itcptmx;itcptmx=tempc;
    tempc=Rmn;Rmn=Rmx;Rmx=tempc;
    for(i=0;i<5;i++){
printf("before:  slpmn: %f  slpmx: %f  slpi: %f\n",slpmn,slpmx,slpp[i]);
       if(slpmn>slpp[i])slpmn=slpp[i];
       if(slpmx<slpp[i])slpmx=slpp[i];
printf("after:  slpmn: %f  slpmx: %f  slpi: %f\n",slpmn,slpmx,slpp[i]);
       if(itcptmn>itcptt[i])itcptmn=itcptt[i];
       if(itcptmx<itcptt[i])itcptmx=itcptt[i];
       if(Rmn>Rr[i])Rmn=Rr[i];
       if(Rmx<Rr[i])Rmx=Rr[i];
      }
    if((slpmn-slpmx)*(slpmn-slpmx)<slpmn*slpmn/800)break;
   }
//  printf("steps: %d\n",step);
//  printf("slope: %f  intercept: %f  residue: %f\n",slpp[0],itcptt[0],rsd[0]);
  printf("%f %f %f %f\n",slpp[0],itcptt[0],Rr[0],rsd[0]);
  return 1;
}
