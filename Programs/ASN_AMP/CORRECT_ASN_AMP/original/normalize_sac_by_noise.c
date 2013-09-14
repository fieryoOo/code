#define MAIN
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "/home/tianye/code/Programs/head/mysac64.h"
#include "/home/tianye/code/Programs/head/sac_db64.h"
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
//  system("pwd\n");
//  if (access(fname,F_OK)==0) {
//     fprintf (stderr,"could be read!!\n");
//     }
//  else {
//     fprintf (stderr," %s could not be read!!!\n",fname);
//     }

  if((fsac = fopen(fname, "rb")) == NULL) {
    printf("could not open sac file to read: %s \n", fname);
    fclose(fsac);
    return 0;
//    abort();
//    exit(1);
  }

  if ( !fsac ) {
    /*fprintf(stderr,"file %s not find\n", fname);*/
    fclose(fsac);
    return NULL;
  }

//  if ( !SHD ) SHD = &SAC_HEADER;

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
SAC_HD tysac;
/*--------------------------------------------------------------*/
float norm_sac(char *fname )
/*--------------------------------------------------------------*/
{
  double dist,minV,maxV,minT,maxT,noiserms2;
  double maxP=40.0;
  double b,e,dt;
  int i,ii,negn,posn,leg,nsamples,n,nn,npow;
  char *fname1, *fltname, *tstr;
  fname1=(char *) malloc(300*sizeof(char));
  fltname = (char *) malloc(300*sizeof(char));
  tstr = (char *) malloc(300*sizeof(char));

  printf("fname is %s\n", fname);

  /*---------------- reading sac file  -------------------*/
  if ( read_sac (fname, sig0, &tysac, SLEN) == NULL ) {
    fprintf(stderr,"file %s not found\n", fname);
    	
    free(fname1);
    free(fltname);
    free(tstr);


    return 0;
  }

  b=tysac.b;
  e=tysac.e;
  if(e-b<2000 || (e<2000 && b>-2000)) 
     {printf("Time length not enough!\n"); return 0;}
  nsamples = tysac.npts;
  dt = tysac.delta;
  negn=(int)(-b/dt);
  leg=negn;
  posn=(int)(e/dt);
  if(posn>leg)leg=posn;
  if( posn+negn+1 != nsamples ){printf("Wrong header file!\n"); return 0;}

n=(int)(1999/dt);
nn=leg-n;
for(i=0;i<=n;i++)
  if(negn+nn+i>=nsamples) sig1[i]=sig0[negn-nn-i];
  else if(negn-nn-i<0) sig1[i]=sig0[negn+nn+i];
  else sig1[i]=(sig0[negn+nn+i]+sig0[negn-nn-i])/2;

  npow=1;

  minV = 2.0;
  maxV = 5.0;
  minT = dist/maxV-maxP/2;
  maxT = dist/minV+maxP;
  if(minT<0)
    minT=0;
  if(maxT>e)
    maxT=e;

double ff1,ff2,ff3,ff4;

  ff1=1./6./1.25;
  ff2=1./6.;
  ff3=1./4.;
  ff4=1./4.*1.25;
  filter4_(&ff1,&ff2,&ff3,&ff4,&npow,&dt,&n,sig1,seis_out);
  noiserms2=0.;
  ii=0;
  int noise_b=(int)(500/dt);
  if((maxT+500)/dt>nn)noise_b=(int)((maxT+500)/dt);
    for(i=noise_b;i<=n-(int)(500/dt);i++) {
//      noise[k]=seis_out[i]*seis_out[i]+noise[k];
      noiserms2=seis_out[i]*seis_out[i]+noiserms2;
      ii++;
    }
    noiserms2=noiserms2/(ii-1);

  for(i=0;i<nsamples;i++) sig0[i]=sig0[i]/noiserms2;
  sprintf(fname1,"%s_norm",fname);
  write_sac(fname1, sig0, &tysac);
  free(fname1);
  free(fltname);
  free(tstr);

  return (noiserms2);
}


/*--------------------------------------------------------------*/
int main (int argc, char *argv[])
/*--------------------------------------------------------------*/
{
  
  FILE *ff;
  int i,nsac;
  char filename[1000000][100], fname[100];
  float noiserms2[1000000];
//  char filename1[300];
  if(argc != 2) 
    {
      printf("Usage: amp_length_normalize list_file\n");
      exit(-1);
    }
  if((ff = fopen(argv[1], "r"))==NULL) {
    printf("Cannot open list file %s to read!\n",argv[1]);
    exit(1);
  }
  for(i=0;;i++) {
    if(fscanf(ff,"%s",&(filename[i][0]))==EOF)
      break;
//    sprintf(filename1,"/utera/tianye/data_ASN_TA/%s",filename);
    noiserms2[i]=norm_sac(filename[i]);
  }
  nsac=i;
  fclose(ff);

  sprintf(fname,"file_noiserms2.lst");
  if((ff = fopen(fname, "w"))==NULL) {
    printf("Cannot open file_noiserms2.lst to write\n");
    exit(1);
  }
  for(i=0;i<nsac;i++) fprintf(ff,"%s  %g\n",filename[i],noiserms2[i]);
  fclose(ff);

return 1;
}
