#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <unistd.h>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;

int main (int argc, char *argv[])
{
  FILE *f;
  char buff[100];
  int i,j,ii,jj,idata,nstep,step,iJ0=3000;
  double datax[5000],datay[5000],rsd[100],dataxx,datayy,rsdmn,rsdfit;
  double J0,J0x[iJ0],J0y[iJ0],A0,A0mx,A0mn,A0stp,A0fit,A0g[100],alpha,alphamx,alphamn,alphastp,alphafit,alphag[100];
  float phvel,omega;

  if(argc != 4)
    {
      printf("Usage: least_squares_line.C [input file] [per] [phvel]\n");
      exit(-1);
    }

  if((f = fopen("/home/tianye/code/Programs/head/J0.txt","r"))==NULL){
     fprintf(stderr,"cannot open J0 file \n");
     return 0;
    }
  for (i=0;i<iJ0;i++) {
     if (fscanf(f,"%lf %lf", &J0x[i],&J0y[i]) != 2) {cout<<"Wrong J0 file<<endl"; return 0;}
    }
  fclose(f);
  phvel=atof(argv[3]);
  omega=2*3.14159265/atof(argv[2]);
  for (i=0;i<iJ0;i++) J0x[i]=J0x[i]*phvel/omega;

  if((f = fopen(argv[1],"r"))==NULL){
     fprintf(stderr,"cannot open data file %s\n",argv[1]);
     return 0;
    }
  for (i=0;;i++) {
     if (fgets(buff,100,f)==NULL)break;
     strtok(buff,"\n");
     if(sscanf(buff,"%lf %lf",&datax[i],&datay[i]) != 2){
        cout<<"Wrong input file!"<<endl; return 0;
       }
//     cout<<datax[i]<<"  "<<datay[i]<<endl;
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

//  A0mx=1.9e-05;A0mn=1.8e-05;
  A0mx=1.0e-4;A0mn=0.;
  alphamx=1.0e-2;alphamn=1.0e-4;
  for(step=0;;step++){
//    printf("A0min: %g  A0max: %g\n",A0mn,A0mx);
    printf("alphamin: %g  alphamax: %g\n",alphamn,alphamx);
    A0stp=(A0mx-A0mn)/9.;
    alphastp=(alphamx-alphamn)/9.;
  
    rsdmn=999999999.;
    for(i=0;i<10;i++){
       A0=A0mn+i*A0stp;
       for(j=0;j<10;j++){
          alpha=alphamn+j*alphastp;
          rsd[i*10+j]=0;
          A0g[i*10+j]=A0;
          alphag[i*10+j]=alpha;
          for(ii=0;ii<idata;ii++){
             for(jj=0;jj<iJ0;jj++)
                if(J0x[jj]>=datax[ii]){J0=J0y[jj];break;}
//             cout<<J0x[jj]<<"  "<<J0<<endl;
             rsd[i*10+j]+=pow((J0*(2*pow(A0,2)/alpha)*exp(-alpha*datax[ii])-datay[ii]),2);
            }
          rsd[i*10+j]=sqrt(rsd[i*10+j]/(idata-1));
         }
      }  
    for(i=0;i<5;i++)
       for(j=i;j<100;j++)
          if(rsd[i]>rsd[j]){
             rsdfit=rsd[i];
             A0fit=A0g[i];
             alphafit=alphag[i];
             rsd[i]=rsd[j];
             A0g[i]=A0g[j];
             alphag[i]=alphag[j];
             rsd[j]=rsdfit;
             A0g[j]=A0fit;
             alphag[j]=alphafit;
            }
  
    A0mn=999999999;A0mx=-999999999;
    alphamn=999999999;alphamx=-999999999;
    for(i=0;i<5;i++){
//printf("residue: %f\n",rsd[i]);
       if(A0mn>A0g[i])A0mn=A0g[i];
       if(A0mx<A0g[i])A0mx=A0g[i];
       if(alphamn>alphag[i])alphamn=alphag[i];
       if(alphamx<alphag[i])alphamx=alphag[i];
      }
    if((alphamn-alphamx)*(alphamn-alphamx)<alphamn*alphamn/10000000)break;
   }
  printf("steps: %d\n",step);
//  printf("slope: %f  intercept: %f  residue: %f\n",slpp[0],itcptt[0],rsd[0]);
  printf("A0: %g  alpha: %g  residue: %g\n",A0g[0],alphag[0],rsd[0]);
  return 1;
}
