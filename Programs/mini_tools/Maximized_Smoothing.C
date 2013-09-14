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

/* CHECK INPUT ARGUMENTS */
  if ( argc < 3 ) {
    fprintf(stderr,"USAGE: Maximize_Smothing [input_sac_file] [maximize_window_#]\n");
    exit(1);
  }

  int i, j, b, e;
  char out_name[100];
  float sig[2000000], sigmax[atoi(arg[2])];
  SAC_HD shd;

  if ( !read_sac (arg[1], sig, &shd, 2000000 ) ) {
    fprintf(stderr,"Cannot read sac file %s\n", arg[1] );
    return 0;
  }

  for(i=0;i<atoi(arg[2]);i++){
     sigmax[i]=0;
     b = (int)floor(i*shd.npts/atof(arg[2])+0.5);
     e = (int)floor((i+1)*shd.npts/atof(arg[2])+0.5);
     for(j=b;j<e;j++)
        if(sig[j]>sigmax[i])sigmax[i]=sig[j];
     for(j=b;j<e;j++)
        sig[j]=sigmax[i];
    }

  sprintf(out_name,"%s.msm\0",arg[1]);
  write_sac (out_name, &(sig[0]), &shd );
}
