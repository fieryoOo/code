#define MAIN
#include "./mysac64.h"
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


SAC_HD shd_cor1_am,shd_cor2_am,shd_cor_temp;
SAC_HD shd_cor1_ph,shd_cor2_ph;

int main(int na, char *arg[])
{ 
  FILE *f1,*f2,*f3;
  int N=2000;
  if(na!=5)
    {
      cout<<"usage:lfstack_from_iris BDH_sac.am BDH_sac.ph BHZ_sac.am BHZ_sac.ph)"<<endl;
      return 0;
    }
  char name[N][6];
  char file_name[300],channel_name[10],year_name[10],month_name[10],month_file_name[300],year_file_name[300];
  float sig1_am[10000000],sig2_am[10000000],sig_temp[10000000];
  float sig1_ph[10000000],sig2_ph[10000000];
  int i,j,k,ii,jj,iii,kk,year,month,path_i;
  double pi=4.0*atan(1.0);

  if(!read_sac (arg[1], sig1_am, &shd_cor1_am, 10000000 ))
    {
      cout<<arg[1]<<"not exist!!"<<endl;
      return 0;
    }
  if(!read_sac (arg[2], sig1_ph, &shd_cor1_ph, 10000000 ))
    {
      cout<<arg[2]<<"not exist!!"<<endl;
      return 0;
    }
  if(!read_sac (arg[3], sig2_am, &shd_cor2_am, 10000000 ))
    {
      cout<<arg[3]<<"not exist!!"<<endl;
      return 0;
    }
  if(!read_sac (arg[4], sig2_ph, &shd_cor2_ph, 10000000 ))
    {
      cout<<arg[4]<<"not exist!!"<<endl;
      return 0;
    }
  int nn=shd_cor1_am.npts;
  if(nn!=shd_cor1_ph.npts)
    {
      fprintf(stderr,"format error!\n");
      return 0;
    }
  if(nn!=shd_cor2_am.npts)
     {
      fprintf(stderr,"format error!\n");
      return 0;
     }
  if(nn!=shd_cor2_ph.npts)
    {
      fprintf(stderr,"format error!\n");
      return 0;
    }
  double dom=shd_cor1_am.delta;
  double b=shd_cor1_am.b;
  double xx,xy,yy,x,y,xy_r,xy_i,x_r,y_r,x_i,y_i;
  double temp_ph,temp_am;
  for(i=50;i<nn-50;i++)
    {
      xx=0;
      //      xy=0;
      yy=0;
      xy_r=0;
      xy_i=0;
      jj=0;
      x_r=0;
      x_i=0;
      y_r=0;
      y_i=0;
      for(ii=-100;ii<=100;ii++)
	{
	  j=ii+i;
	  x=sig1_am[j];
	  y=sig2_am[j];
	  xx+=x*x;
	  yy+=y*y;
	  // xy+=x*y;//*cos(sig1_ph[j]-sig2_ph[j]);
	  xy_r+=x*y*cos(sig1_ph[j]-sig2_ph[j]);
	  xy_i+=x*y*sin(sig1_ph[j]-sig2_ph[j]);
	  jj++;
	}
      temp_ph=atan(xy_i/xy_r);
      if(xy_r<0)
	temp_ph+=pi;
      if(temp_ph>pi)
	temp_ph-=2*pi;
      temp_am=sqrt(xy_r*xy_r+xy_i*xy_i);
      cout<<dom*i+b<<" "<<temp_am/sqrt(xx*yy)<<" "<<temp_am/xx<<" "<<temp_ph<<endl;
      //" "<<sqrt(xx)<<" "<<sqrt(yy)<<" "<<temp_am<<endl;
      //cout<<dom*i+b<<" "<<sig2_am[i]/sig1_am[i]<<endl;
      if(dom*i+b>1)
	break;
    }
}
