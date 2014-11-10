#define MAIN
#include "../inc/mysac64.h"
#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace std;


/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
        SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
/*..........................................................................*/
        if((fsac=fopen(fname, "rb")) == NULL) return NULL;

        if ( !SHD ) SHD = &SAC_HEADER;

         fread(SHD,sizeof(SAC_HD),1,fsac);

         if ( SHD->npts > nmax ) {
           fprintf(stderr,
           "ATTENTION !!! in the file %s npts exceeds limit  %d",fname,nmax);
           SHD->npts = nmax;
         }

         fread(sig,sizeof(float),(int)(SHD->npts),fsac);

         fclose (fsac);

   /*-------------  calculate from t0  ----------------*/
   {
        int eh, em ,i;
        float fes;
        char koo[9];

        for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
        koo[8] = '\0';

        SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
         SHD->nzsec + SHD->nzmsec*.001;

        sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

        SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}

        return SHD;
}




/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
        void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
 int i;
/*..........................................................................*/
        fsac = fopen(fname, "wb");

        if ( !SHD ) SHD = &SAC_HEADER;


        SHD->iftype = (int)ITIME;
        SHD->leven = (int)TRUE;

        SHD->lovrok = (int)TRUE;
        SHD->internal4 = 6L;



  /*+++++++++++++++++++++++++++++++++++++++++*/
     SHD->depmin = sig[0];
     SHD->depmax = sig[0];
 
   for ( i = 0; i < SHD->npts ; i++ ) {
    if ( SHD->depmin > sig[i] ) SHD->depmin = sig[i];
    if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }
   fwrite(SHD,sizeof(SAC_HD),1,fsac);

   fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);


   fclose (fsac);
}


SAC_HD shd_cor_EE;
SAC_HD shd_cor_EN;
SAC_HD shd_cor_NE;
SAC_HD shd_cor_NN;
SAC_HD shd_cor_ZE;
SAC_HD shd_cor_ZN;
SAC_HD shd_cor_NZ;
SAC_HD shd_cor_EZ;
SAC_HD shd_cor_temp;

int main(int na, char *arg[])
{ 
  FILE *f1,*f2,*f3;
  int N=5000;
  double PI;  
  PI=atan(1.0)*4;
  double temp1,temp2;
  double cos1,cos2,sin1,sin2;
  char name[N][6];
  char name_ZZ[200],name_ZE[200],name_ZN[200];
  char name_EZ[200],name_EE[200],name_EN[200];
  char name_NZ[200],name_NE[200],name_NN[200];
  char name_ZR[200],name_ZT[200];
  char name_RZ[200],name_RR[200],name_RT[200];
  char name_TZ[200],name_TR[200],name_TT[200];
  float sig_ZZ[100000],sig_ZE[100000],sig_ZN[100000];
  float sig_EZ[100000],sig_EE[100000],sig_EN[100000];
  float sig_NZ[100000],sig_NE[100000],sig_NN[100000];
  float sig_ZR[100000],sig_ZT[100000];
  float sig_RZ[100000],sig_RR[100000],sig_RT[100000];
  float sig_TZ[100000],sig_TR[100000],sig_TT[100000]; 
  float sig_temp[100000];
  if(na!=5)
    {
      cout<<"usage:roate_ZNE2ZRT stalist path ns1 ns2"<<endl;
      return 0;
    }
  f1=fopen(arg[1],"r");
  int i,j,k,ii,jj,ns1,ns2;
  for(i=0;i<N;i++)
    {
      if(fscanf(f1,"%s",name[i])==EOF) break;
      cout<<name[i]<<endl;
    }
  fclose(f1);

  cout<<"number of stations read "<<i<<endl;
  char buff[300];
  //  sprintf(buff,"if [ ! -d STACK_ZRT ]; then mkdir STACK_ZRT; fi");
  //system(buff);
  //  sprintf(buff,"if [ ! -d STACK_TT ]; then mkdir STACK_TT; fi");
  // system(buff);
  //sprintf(buff,"if [ ! -d STACK_RR ]; then mkdir STACK_RR; fi");
  //system(buff);
  //sprintf(buff,"if [ ! -d STACK_TR ]; then mkdir STACK_TR; fi");
  //system(buff);
  //sprintf(buff,"if [ ! -d STACK_RT ]; then mkdir STACK_RT; fi");
  //system(buff);

  // section the outer loop
  ns1=atoi(arg[3]);
  ns2=atoi(arg[4]);
  if(ns1==0 && ns2==0)
    {
      ns1=1;
      ns2=i;
    }
  if(ns1<1)
    {
      ns1=1;
    }
  if(ns2>i)
    {
      ns2=i;
    }
  ns1=ns1-1;
  ns2=ns2-1;
  printf("%d %d\n",ns1,ns2);


  for(j=ns1;j<=ns2;j++)
  //for(j=0;j<i-1;j++)   Fan-Chi for H/V
    {
      //sprintf(buff,"if [ ! -d STACK_ZRT/%s ]; then mkdir STACK_ZRT/%s; fi",name[j],name[j]);
      //system(buff);
      //      sprintf(buff,"if [ ! -d STACK_TT/%s ]; then mkdir STACK_TT/%s; fi",name[j],name[j]);
      //system(buff);
      //sprintf(buff,"if [ ! -d STACK_RR/%s ]; then mkdir STACK_RR/%s; fi",name[j],name[j]);
      //system(buff);
      //sprintf(buff,"if [ ! -d STACK_TR/%s ]; then mkdir STACK_TR/%s; fi",name[j],name[j]);
      //system(buff);
      //sprintf(buff,"if [ ! -d STACK_RT/%s ]; then mkdir STACK_RT/%s; fi",name[j],name[j]);
      //system(buff);

      for(k=0;k<i;k++)
      //for(k=j+1;k<i;k++)  Fan-Chi for H/V
	{
	  cout<<"working on "<<name[j]<<" "<<name[k]<<endl;
	  sprintf(name_ZZ,"%s/%s/COR_%s_%s.SAC_ZZ",arg[2],name[j],name[j],name[k]);
	  sprintf(name_ZE,"%s/%s/COR_%s_%s.SAC_ZE",arg[2],name[j],name[j],name[k]);
	  sprintf(name_ZN,"%s/%s/COR_%s_%s.SAC_ZN",arg[2],name[j],name[j],name[k]);
	  sprintf(name_EZ,"%s/%s/COR_%s_%s.SAC_EZ",arg[2],name[j],name[j],name[k]);
	  sprintf(name_EE,"%s/%s/COR_%s_%s.SAC_EE",arg[2],name[j],name[j],name[k]);
	  sprintf(name_EN,"%s/%s/COR_%s_%s.SAC_EN",arg[2],name[j],name[j],name[k]);
	  sprintf(name_NZ,"%s/%s/COR_%s_%s.SAC_NZ",arg[2],name[j],name[j],name[k]);
	  sprintf(name_NE,"%s/%s/COR_%s_%s.SAC_NE",arg[2],name[j],name[j],name[k]);
	  sprintf(name_NN,"%s/%s/COR_%s_%s.SAC_NN",arg[2],name[j],name[j],name[k]);
	  sprintf(name_ZR,"%s/%s/COR_%s_%s.SAC_ro_ZR",arg[2],name[j],name[j],name[k]);
	  sprintf(name_ZT,"%s/%s/COR_%s_%s.SAC_ro_ZT",arg[2],name[j],name[j],name[k]);
	  sprintf(name_RZ,"%s/%s/COR_%s_%s.SAC_ro_RZ",arg[2],name[j],name[j],name[k]);
	  sprintf(name_RR,"%s/%s/COR_%s_%s.SAC_ro_RR",arg[2],name[j],name[j],name[k]);
	  sprintf(name_RT,"%s/%s/COR_%s_%s.SAC_ro_RT",arg[2],name[j],name[j],name[k]);
	  sprintf(name_TZ,"%s/%s/COR_%s_%s.SAC_ro_TZ",arg[2],name[j],name[j],name[k]);
	  sprintf(name_TR,"%s/%s/COR_%s_%s.SAC_ro_TR",arg[2],name[j],name[j],name[k]);
	  sprintf(name_TT,"%s/%s/COR_%s_%s.SAC_ro_TT",arg[2],name[j],name[j],name[k]);
	  if(read_sac (name_TT, sig_temp, &shd_cor_temp, 100000 ))
            {
              cout<<name_TT<<" exist!!"<<endl;
//              continue;
            }
	  if(!read_sac (name_EE, sig_temp, &shd_cor_temp, 100000 ))
            {
              cout<<name_EE<<" not exist!!"<<endl;
              continue;
            }
	  

	  f2=fopen("runsac.csh","w");
          fprintf(f2,"sac <<END\n");
          fprintf(f2,"r %s %s %s %s %s %s %s %s %s\n",name_ZZ,name_ZE,name_ZN,name_EZ,name_EE,name_EN,name_NZ,name_NE,name_NN);
	  fprintf(f2,"w %s %s %s %s %s %s %s %s %s\n",name_ZZ,name_ZE,name_ZN,name_EZ,name_EE,name_EN,name_NZ,name_NE,name_NN);
          fprintf(f2,"quit\n");
	  fprintf(f2,"END\n\n");
          fclose(f2);
          system("csh runsac.csh");
	  if(!read_sac (name_ZE, sig_ZE, &shd_cor_ZE, 100000 ))
	    {
	      cout<<name_ZE<<" not found!!"<<endl;
	      continue;
	    }
	  if(!read_sac (name_ZN, sig_ZN, &shd_cor_ZN, 100000 ))
	    {
	      cout<<name_ZN<<" not found!!"<<endl;
	      continue;
	    }
	  if(!read_sac (name_NZ, sig_NZ, &shd_cor_NZ, 100000 ))
	    {
	      cout<<name_NZ<<" not found!!"<<endl;
	      continue;
	    }	
	  if(!read_sac (name_EZ, sig_EZ, &shd_cor_EZ, 100000 ))
	    {
	      cout<<name_EZ<<" not found!!"<<endl;
	      continue;
	    }
	  if(!read_sac (name_EE, sig_EE, &shd_cor_EE, 100000 ))
	    {
	      cout<<name_EE<<" not found!!"<<endl;
	      continue;
	    }
	  if(!read_sac (name_EN, sig_EN, &shd_cor_EN, 100000 ))
	    {
	      cout<<name_EN<<" not found!!"<<endl;
	      continue;
	    }
	  if(!read_sac (name_NE, sig_NE, &shd_cor_NE, 100000 ))
	    {
	      cout<<name_NE<<" not found!!"<<endl;
	      continue;
	    }	
	  if(!read_sac (name_NN, sig_NN, &shd_cor_NN, 100000 ))
	    {
	      cout<<name_NN<<" not found!!"<<endl;
	      continue;
	    }
	  temp1=fabs(cos((shd_cor_NN.cmpaz-shd_cor_EE.cmpaz)/180.0*PI));
	  temp2=fabs(cos((shd_cor_NN.user9-shd_cor_EE.user9)/180.0*PI));
	  if(temp1>0.1||temp2>0.1)
	    {
	      cout<<"cmpaz incompatible!!"<<endl;
	      continue;
	    }
//	  temp1=(shd_cor_NN.az-shd_cor_NN.cmpaz-180)/180.0*PI;
	  temp1=(shd_cor_NN.baz-shd_cor_NN.user9)/180.0*PI;
	  cos1=cos(temp1);
	  sin1=sin(temp1);
//	  temp2=(shd_cor_NN.baz-shd_cor_NN.user9-180)/180.0*PI;
	  temp2=(shd_cor_NN.az-shd_cor_NN.cmpaz)/180.0*PI;
	  cos2=cos(temp2);
	  sin2=sin(temp2);
	  //cout<<temp1<<" "<<sin1<<" "<<cos1<<" "<<temp2<<" "<<sin2<<" "<<cos2<<endl;
	  for(ii=0;ii<shd_cor_NN.npts;ii++)
	    {
//	      sig_ZR[ii]=-1*sin2*sig_ZE[ii]-cos2*sig_ZN[ii];
	      sig_ZR[ii]=1*sin2*sig_ZE[ii]-cos2*sig_ZN[ii];
	      sig_ZT[ii]=-1*cos2*sig_ZE[ii]+sin2*sig_ZN[ii];

//	      sig_RZ[ii]=1*sin1*sig_EZ[ii]+cos1*sig_NZ[ii];
	      sig_RZ[ii]=-1*sin1*sig_EZ[ii]+cos1*sig_NZ[ii];
	      sig_TZ[ii]=1*cos1*sig_EZ[ii]-sin1*sig_NZ[ii];

	      sig_TT[ii]=-1*cos1*cos2*sig_EE[ii]+cos1*sin2*sig_EN[ii]-sin1*sin2*sig_NN[ii]+sin1*cos2*sig_NE[ii];
	      sig_RR[ii]=-1*sin1*sin2*sig_EE[ii]-sin1*cos2*sig_EN[ii]-cos1*cos2*sig_NN[ii]-cos1*sin2*sig_NE[ii];
//	      sig_TR[ii]=-1*cos1*sin2*sig_EE[ii]-cos1*cos2*sig_EN[ii]+sin1*cos2*sig_NN[ii]+sin1*sin2*sig_NE[ii];
	      sig_TR[ii]=-1*cos1*sin2*sig_EE[ii]-cos1*cos2*sig_EN[ii]+sin1*cos2*sig_NN[ii]-sin1*sin2*sig_NE[ii];
              sig_RT[ii]=-1*sin1*cos2*sig_EE[ii]-sin1*sin2*sig_EN[ii]+cos1*sin2*sig_NN[ii]-cos1*cos2*sig_NE[ii];
	    }
	  write_sac (name_TT, sig_TT, &shd_cor_NN);
	  write_sac (name_RR, sig_RR, &shd_cor_NN);
	  write_sac (name_TR, sig_TR, &shd_cor_NN);
	  write_sac (name_RT, sig_RT, &shd_cor_NN);
	  write_sac (name_ZR, sig_ZR, &shd_cor_NN);
	  write_sac (name_ZT, sig_ZT, &shd_cor_NN);
	  write_sac (name_RZ, sig_RZ, &shd_cor_NN);
	  write_sac (name_TZ, sig_TZ, &shd_cor_NN);
	}
	  

    }

 
}
