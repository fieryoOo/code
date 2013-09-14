#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

int sort ( double *datax, double *datay1, double *datay2, double *datay3, double *datay4, int ndat )
{
  int i, j, ii;
  double datamin, dataxx, datayy;

  for(i=0;i<ndat;i++) {
     datamin=datax[i]; ii=i;
     for(j=i+1;j<ndat;j++)
        if(datamin>datax[j]){
           datamin=datax[j];
           ii=j;
        }
     if(ii==i) continue;
     dataxx=datax[i];
     datax[i]=datax[ii];
     datax[ii]=dataxx;
     datayy=datay1[i];
     datay1[i]=datay1[ii];
     datay1[ii]=datayy;
     datayy=datay2[i];
     datay2[i]=datay2[ii];
     datay2[ii]=datayy;
     datayy=datay3[i];
     datay3[i]=datay3[ii];
     datay3[ii]=datayy;
     datayy=datay4[i];
     datay4[i]=datay4[ii];
     datay4[ii]=datayy;
     //datayy=datay5[i];
     //datay5[i]=datay5[ii];
     //datay5[ii]=datayy;
  }

  return 1;
}

int sort2 ( double *datax, double *datay1, int ndat )
{
  int i, j, ii;
  double datamin, dataxx, datayy;

  for(i=0;i<ndat;i++) {
     datamin=datax[i]; ii=i;
     for(j=i+1;j<ndat;j++)
        if(datamin>datax[j]){
           datamin=datax[j];
           ii=j;
        }
     if(ii==i) continue;
     dataxx=datax[i];
     datax[i]=datax[ii];
     datax[ii]=dataxx;
     datayy=datay1[i];
     datay1[i]=datay1[ii];
     datay1[ii]=datayy;
  }

  return 1;
}

int fit_line ( double *dataxo, double *datayo, int idatao, double *slpout, double *itcptout, double *sdout )
{
  FILE *f;
  int i,j,ii,nstep,step, idata;
  double slp,slpmx,slpmn,slpstp,itcpt,itcptmx,itcptmn,itcptstp,slpfit,itcptfit,slpp[100],itcptt[100],rsd[100],dataxx,datayy,datamin,rsdmn,rsdfit;
  double datax[idatao], datay[idatao];


  for(i=0;i<idatao;i++) {
     datamin=dataxo[i]; ii=i;
     for(j=i+1;j<idatao;j++)
        if(datamin>dataxo[j]){
           datamin=dataxo[j];
           ii=j;
        }
     if(ii==i) continue;
     dataxx=dataxo[i];
     datayy=datayo[i];
     dataxo[i]=dataxo[ii];
     datayo[i]=datayo[ii];
     dataxo[ii]=dataxx;
     datayo[ii]=datayy;
  }
  
  idata=0;
  for(i=0;i<idatao;i++) {
     if(datayo[i]==-1) continue;
     datax[idata] = dataxo[i];
     datay[idata] = datayo[i];
     idata++;
  }

  slpmx=-999999999.;slpmn=999999999.;
  itcptmx=-999999999.;itcptmn=999999999.;
  dataxx=(datax[0]+datax[1]+datax[2])/3;
  datayy=(datay[0]+datay[1]+datay[2])/3;
  for(i=idata/2;i<idata;i++){
     slp=(datay[i]-datayy)/(datax[i]-dataxx);
     //if(slp>0) continue;
     if(slpmx<slp)slpmx=slp;
     if(slpmn>slp)slpmn=slp;
    }
  dataxx=(datax[idata-1]+datax[idata-2]+datax[idata-3])/3;
  datayy=(datay[idata-1]+datay[idata-2]+datay[idata-3])/3;
  for(i=0;i<idata/2;i++){
     slp=(datayy-datay[i])/(dataxx-datax[i]);
     //if(slp>0) continue;
     if(slpmx<slp)slpmx=slp;
     if(slpmn>slp)slpmn=slp;
    }
  for(i=0;i<idata;i++){
     itcpt=datay[i]-datax[i]*slpmx;
     if(itcptmn>itcpt)itcptmn=itcpt;
     itcpt=datay[i]-datax[i]*slpmn;
     if(itcptmx<itcpt)itcptmx=itcpt;
    }

 // nstep=atoi(argv[3]);

  for(step=0;;step++){
//    printf("slpmin: %f  slpmax: %f\n",slpmn,slpmx);
    slpstp=(slpmx-slpmn)/9;
    itcptstp=(itcptmx-itcptmn)/9;

    rsdmn=999999999.;
    for(i=0;i<10;i++){
       slp=slpmn+i*slpstp;
       for(j=0;j<10;j++){
          itcpt=itcptmn+j*itcptstp;
          rsd[i*10+j]=0;
          slpp[i*10+j]=slp;
          itcptt[i*10+j]=itcpt;
          for(ii=0;ii<idata;ii++){
             rsd[i*10+j]+=(slp*datax[ii]+itcpt-datay[ii])*(slp*datax[ii]+itcpt-datay[ii]);
            }
          rsd[i*10+j]=sqrt(rsd[i*10+j]/(idata-1));
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
   }
//  printf("steps: %d\n",step);
//  printf("slope: %f  intercept: %f  residue: %f\n",slpp[0],itcptt[0],rsd[0]);
//  printf("%f %f %f\n",slpp[0],itcptt[0],rsd[0]);
  *slpout=slpp[0]; *itcptout=itcptt[0]; *sdout=rsd[0];
  return 1;
}
