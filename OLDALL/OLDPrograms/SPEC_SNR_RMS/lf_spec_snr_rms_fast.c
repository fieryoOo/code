#define MAIN
#include "/home/tianye/code/Programs/head/mysac64.h"
#include "/home/tianye/code/Programs/head/sac_db64.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "/home/tianye/code/Programs/head/koftan.h"
#include "/home/tianye/code/Programs/head/gl_const.h"
#include "/home/tianye/code/Programs/head/mymacro.h"

#define SLEN 200000

void filter4_(double *f1,double *f2,double *f3,double *f4,int *npow,
              double *dt,int *n, float seis_in[], float seis_out[]);


/*--------------------------------------------------------------------------*/
SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*--------------------------------------------------------------------------*/
/* function to read sac files given the name, fname. The function outputs the time signal to the pointer sig, fills the header SHD, if the signal has fewer than nmax points */
{
  FILE *fsac;

  if((fsac = fopen(fname, "rb")) == NULL) {
    printf("could not open sac file to read%s \n", fname);
    exit(1);
  }

  if ( !fsac ) {
    /*fprintf(stderr,"file %s not find\n", fname);*/
    return NULL;
  }

  if ( !SHD ) SHD = &SAC_HEADER;

  fread(SHD,sizeof(SAC_HD),1,fsac);

  if ( SHD->npts > nmax ) {
    fprintf(stderr,"ATTENTION !!! %s npts is limited to %d.\n",fname,nmax);
    SHD->npts = nmax;
  }

  fread(sig,sizeof(float),(int)(SHD->npts),fsac);
  fclose (fsac);

/*-------------  calcule de t0  ----------------*/
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


/*--------------------------------------------------------------------------*/
void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*--------------------------------------------------------------------------*/
{
  FILE *fsac;
  int i;
  if((fsac = fopen(fname, "wb"))==NULL) {
    printf("could not open sac file to write\n");
    exit(1);
  }

  if ( !SHD ) {
    SHD = &SAC_HEADER;
  }

  SHD->iftype = (int)ITIME;
  SHD->leven = (int)TRUE;
  SHD->lovrok = (int)TRUE;
  SHD->internal4 = 6L;
  SHD->depmin = sig[0];
  SHD->depmax = sig[0];

  for ( i = 0; i < SHD->npts ; i++ ) {
    if ( SHD->depmin > sig[i] ) {
      SHD->depmin = sig[i];
    }
    if ( SHD->depmax < sig[i] ) {
      SHD->depmax = sig[i];
    }
   }

  fwrite(SHD,sizeof(SAC_HD),1,fsac);
  fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);

  fclose (fsac);
}


float sig0[SLEN],seis_out[SLEN];
char fname[300];
/*--------------------------------------------------------------*/
int get_snr(char *fname )
/*--------------------------------------------------------------*/
{
  FILE *fp1, *f2;
  double dist,minV,maxV,minT,maxT,window,signalmax,noisemax,noiserms;
  double maxP=150.0;
  int nf=32;
  double per[nf],f[nf],snr[nf],amp_max[nf];
  double minP=5.0, num;
  double b,e,fb,fe,step;
  double per1, per2, per3, per4;
  int k,i, j,ii;
  char fname1[300], fltname[200];

  printf("fname is %s\n", fname);

  /*---------------- reading sac file  -------------------*/
  if ( read_sac (fname, sig0, &SAC_HEADER, SLEN) == NULL ) {
    fprintf(stderr,"file %s not found\n", fname);
    return 0;
  }

  fb=1.0/maxP;
  fe=1.0/minP;
  step=(log(fb)-log(fe))/(nf-1);
  for(k=0;k<nf;k++) {
    f[k]=exp(log(fe)+k*step);
    per[k]=1/f[k];
  }
  b=SAC_HEADER.b;
  e=SAC_HEADER.e;
  dist = SAC_HEADER.dist;
  
  double dt;
  int nsamples,n,npow;
  
  dt = SAC_HEADER.delta;
  nsamples = SAC_HEADER.npts;

  npow=1;
  n=nsamples;

  minV = 2.0;
  maxV = 5.0;
  minT = dist/maxV-maxP;
  maxT = dist/minV+2*maxP;
  if(minT<0)
    minT=0;
  if(maxT>e)
    maxT=e;
  window=maxT-minT;

  double ff1,ff2,ff3,ff4;
  char buff[300];
  int imk;
  for(k=1;k<nf-1;k++) 
    {
      if((f[k+1] != 0) && (f[k-1] != 0))
	{
	  ff1 = f[k+1]/1.25;
	  ff2 = f[k+1];
	  ff3 = f[k-1];
	  ff4 = f[k-1]*1.25;
	}
      else return 1;
      filter4_(&ff1,&ff2,&ff3,&ff4,&npow,&dt,&n,sig0,seis_out);
      // if(ff2<0.01&&ff3>0.01)
      //{
      //  sprintf(buff,"test.sac");
      //  write_sac(buff, seis_out, &SAC_HEADER);
      //}
      
      //fprintf(stderr,"%lf %lf %lf %lf\n", ff1, ff2, ff3, ff4);

    ///////////////////////////////////////////////////////////
    // bp filtered signal has been obtained, now go through and 
    // find the snr
    ///////////////////////////////////////////////////////////
    signalmax=0;
    noisemax=0;
    //for(i=int(minT);i<maxT;i++)
    for(i=(int)minT-b;i<maxT-b;i++) {
      if(seis_out[i] < 0) num = seis_out[i]*(-1);
      else num = seis_out[i];
      if(num>signalmax) { 
        imk=i; signalmax=num; }
    }
    //fprintf(stderr,"%lf %lf %lf\n", minT, maxT, signalmax);
    noiserms=0;
    ii=0;
    //printf("imarker: %d\n",imk);
    for(i=imk+400;i<e-b-100;i++) {
      noiserms=seis_out[i]*seis_out[i]+noiserms;
      ii++;
      //if(seis_out[i] < 0) num = seis_out[i]*(-1);
      //else num = seis_out[i];
      //if(num>noisemax)
      //noisemax=num;
    }
    noiserms=sqrt(noiserms/(ii-1));
  //  printf("%f %f %f %lf %lf\n", 1/f[k], ff2, ff3, signalmax, noiserms);
    //fprintf(stderr,"%lf\n",noiserms);
    snr[k]=signalmax/noiserms;
    amp_max[k]=signalmax;
    //    snr[k]=signalmax/noisemax;      
  }
  sprintf(fname1,"%s_snr_rms.txt",fname);
  f2=fopen(fname1,"w");
  for(i=1;i<nf-1;i++) {
    fprintf(f2,"%lf %lf %lf\n",per[i],snr[i],amp_max[i]);
  }
  printf("maxP=%f minP=%f maxV=%f minV=%f!!!\n",maxP,minP,maxV,minV);
  fclose(f2);
  return 1;
}


/*--------------------------------------------------------------*/
int main (int argc, char *argv[])
/*--------------------------------------------------------------*/
{
  
  FILE *ff;
  char filename[200];
  if(argc != 2) 
    {
      printf("Usage: spec_snr_rms list_file\n");
      exit(-1);
    }
  if((ff = fopen(argv[1], "r"))==NULL) {
    printf("Cannot open list.lst file. A file containing the names of");
    printf(" all cut files for getting spectral SNR is required to run spectral_snr*.\n");
    exit(1);
  }

  for(;;) {
    if(fscanf(ff,"%s",&filename)==EOF)
      break;
    get_snr(filename);
  }
}
