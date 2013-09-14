#define MAIN
#include "/home/tianye/code/Programs/head/mysac64.h"
#include "/home/tianye/code/Programs/head/sac_db64.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SLEN 200000

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

int reverse(char *fname,char *cen_sta)
{
FILE *f1;
int N=3000;
int i,npts;
double b,e,stla,stlo,evla,evlo;
float sig0[SLEN],sig1[SLEN];
SAC_HD tysac;
char evnm[16],stnm[16];
char out_name[300];

if( read_sac (fname, sig0, &tysac, SLEN ) == NULL) {
    fprintf(stderr,"could not open sac file to read: %s \n", fname);
    return 0;
  }

b=tysac.b;
e=tysac.e;
strcpy(evnm,tysac.kevnm);
strcpy(stnm,tysac.kstnm);
stla=tysac.stla;
stlo=tysac.stlo;
evla=tysac.evla;
evlo=tysac.evlo;
npts=tysac.npts;
for(i=0;i<npts;i++) {
   sig1[i]=sig0[npts-1-i];
  }
tysac.b=-e;
tysac.e=-b;
tysac.stla=evla;
tysac.stlo=evlo;
tysac.evla=stla;
tysac.evlo=stlo;
strcpy(tysac.kevnm,stnm);
strcpy(tysac.kstnm,evnm);

sprintf(out_name,"COR_%s_%s.SAC",stnm,evnm);
write_sac(out_name,sig1,&tysac);

return 1;
}

int main (int argc, char *argv[])
{
FILE *ff;
char filename[300],sta[10];
if(argc != 3)
  {
   printf("Usage: sort_COR_by_centre_sta.C [centre_sta] [input_SAC_lst]\n");
   exit(-1);
  }
if((ff = fopen(argv[2], "r"))==NULL) {
   printf("Cannot open input sac file list %s.\n",argv[2]);
   exit(1);
  }

for(;;) {
   if(fscanf(ff,"%s",&(filename[0]))==EOF)
   break;
   sscanf(filename,"COR_%[^'_']",&sta);
   if(!strcmp(sta,argv[1]))continue;
   reverse(filename,argv[1]);
   printf("%s reversed...\n",filename);
  }
fclose(ff);
return 1;
}
