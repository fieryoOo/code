#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "/home/tianye/code/Programs/head/mysac64.h"
#include "/home/tianye/code/Programs/head/64_sac_db.h"
#include <string.h>

SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax);
void write_sac (char *fname, float *sig, SAC_HD *SHD);

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
    /*fprintf(stderr,"file %s not found\n", fname);*/
    return NULL;
  }

  if ( !SHD ) SHD = &SAC_HEADER;

  fread(SHD,sizeof(SAC_HD),1,fsac);

  if ( SHD->npts > nmax ) {
    fprintf(stderr,"ATTENTION !!! dans le fichier %s npts est limite a %d",fname,nmax);
    SHD->npts = nmax;
  }

  fread(sig,sizeof(float),(int)(SHD->npts),fsac);
  fclose (fsac);

/*-------------  calcule de t0  ----------------*/
   {
     int eh, em ,i;
     float fes;
     char koo[9];

     for ( i = 0; i < 8; i++ ) {
       koo[i] = SHD->ko[i];
     }
     koo[8] = 0;

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

/*--------------------------------------------------------------------------*/
int main (int argc, char *arg[])
/*--------------------------------------------------------------------------*/
{
  int i, j, b, e, window_b, window_e, h_l;
  char out_name[100];
  float sig[2000000], window_sum[2000000];
  SAC_HD shd;

/* CHECK INPUT ARGUMENTS */
  if ( argc < 5 ) {
    fprintf(stderr,"USAGE: Whitening [input_sac_file] [lower_end] [upper_end] [half_window_n]\n");
    return 0;
  }

  if ( !read_sac (arg[1], sig, &shd, 2000000 ) ) {
    fprintf(stderr,"Cannot read sac file %s\n", arg[1] );
    return 0;
  }

  b=(int)ceil(atof(arg[2])/shd.delta);
  e=(int)ceil(atof(arg[3])/shd.delta);
  h_l=atoi(arg[4]);
  for(i=0;i<b;i++)sig[i]=0;
  for(i=e;i<shd.npts;i++)sig[i]=0;

  for(i=b;i<e;i++){
      window_sum[i]=0;
      window_b=i-h_l;
      if(window_b<b)window_b=b;
      window_e=i+h_l;
      if(window_e>e)window_e=e;
      for(j=window_b;j<window_e;j++)window_sum[i]+=fabs(sig[j]);
      window_sum[i]/=(window_e-window_b);
     }
  for(i=b;i<e;i++) sig[i]=sig[i]/window_sum[i];

  sprintf(out_name,"%s.wt\0",arg[1]);
  write_sac (out_name, &(sig[0]), &shd );
}
