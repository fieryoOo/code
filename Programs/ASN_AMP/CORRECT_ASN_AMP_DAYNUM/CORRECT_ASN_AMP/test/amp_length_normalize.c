#define MAIN
#include "/home/linf/NOISE_CODA_64/mysac64.h"
#include "/home/linf/NOISE_CODA_64/sac_db64.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "/home/linf/PROGS/HEADER_FILES/koftan.h"
#include "/home/linf/PROGS/HEADER_FILES/gl_const.h"
#include "/home/linf/PROGS/HEADER_FILES/mymacro.h"

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


float sig0[SLEN],seis_out[SLEN],sig1[SLEN];
char fname[300];
/*--------------------------------------------------------------*/
int get_snr(char *fname )
/*--------------------------------------------------------------*/
{
  FILE *fp1, *f2;
  double dist,minV,maxV,minT,maxT,window,signalmax,noisemax;
  double maxP=30.0;
  int nf=8;
  double per[nf],f[nf],snr[nf],amp_max[nf],noiserms[nf];
  double minP=5.0, num;
  double b,e,fb,fe,step;
  double per1, per2, per3, per4;
  int k,i, j,ii,bi,ei;
  char fname1[300], fltname[200], tstr[300];

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
//  e=-b;
  bi=(int)b;
  ei=(int)e;
//printf("b: %f  e: %f\n",b,e);
//printf("bi: %d  ei: %d\n",bi,ei);
for(i=-bi;i<=-bi+ei;i++) {
  
  sig1[i+bi]=(sig0[i]+sig0[-bi*2-i])/2;
//if(sig0[i]*sig0[i]>3)printf("i: %d  i: %f  -i: %f  s: %f\n",i+bi,sig0[i],sig0[-bi*2-i],sig1[i+bi]);
  }
  dist = SAC_HEADER.dist;
//printf("%f\n",dist);  
  double dt;
  int nsamples,n,npow;
  
  dt = SAC_HEADER.delta;
  nsamples = SAC_HEADER.npts;

  npow=1;
  n=nsamples;

  minV = 2.0;
  maxV = 3.5;
  minT = dist/maxV-maxP/2;
  maxT = dist/minV+maxP;
  if(minT<0)
    minT=0;
  if(maxT>e)
    maxT=e;
  window=maxT-minT;
//printf("%f  %f\n",minT,maxT);

  double ff1,ff2,ff3,ff4;
  char buff[300];
  for(k=1;k<nf-1;k++) 
    {
      if((f[k+1] != 0) && (f[k-1] != 0))
	{
	  ff1 = f[k+1]/1.25;
	  ff2 = f[k+1];
	  ff3 = f[k-1];
	  ff4 = f[k-1]*1.25;
//printf("%f %f %f %f\n",1/ff1,1/ff2,1/ff3,1/ff4);
	}
      else return 1;
      filter4_(&ff1,&ff2,&ff3,&ff4,&npow,&dt,&n,sig1,seis_out);
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

    for(i=(int)minT;i<maxT;i++) {
      if(seis_out[i] < 0) num = seis_out[i]*(-1);
      else num = seis_out[i];
      if(num>signalmax) 
      signalmax=num;
    }
//if(signalmax>3)printf("signalmax: %f\n",signalmax);
    //fprintf(stderr,"%lf %lf %lf\n", minT, maxT, signalmax);
    noiserms[k]=0;
    ii=0;
int noise_b=1500;
if(maxT+500>1500)noise_b=(int)maxT+500;
    for(i=noise_b;i<2500;i++) {
      noiserms[k]=seis_out[i]*seis_out[i]+noiserms[k];
      ii++;
      //if(seis_out[i] < 0) num = seis_out[i]*(-1);
      //else num = seis_out[i];
      //if(num>noisemax)
      //noisemax=num;
    }
    noiserms[k]=sqrt(noiserms[k]/(ii-1));
    //fprintf(stderr,"%lf\n",noiserms);
//printf("signalmax: %f  noiserms: %f",signalmax,noiserms);
    snr[k]=signalmax/noiserms[k];
    amp_max[k]=signalmax;
    //    snr[k]=signalmax/noisemax;      
  }
  sprintf(fname1,"%s_test_snr.txt",fname);
  f2=fopen(fname1,"w");
  if(snr[(int)(nf/5)]<15 || snr[nf-1-(int)(nf/5)]<15){
    fprintf(stderr,"SNR of %s less than 15. skipped!\n",fname );
    sprintf(tstr,"rm -f %s\n",fname1);
    system(tstr);
    return 0;
  }
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
      printf("Usage: amp_length_normalize list_file\n");
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
