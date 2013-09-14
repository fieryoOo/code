#define MAIN
#include "/home/linf/NOISE_CODA_64/mysac64.h"
#include "/home/linf/NOISE_CODA_64/sac_db64.h"
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
  if((fsac = fopen(fname, "rb")) == NULL) {
    printf("could not open sac file to read: %s \n", fname);
    fclose(fsac);
    return 0;
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

int main(int argc, char *argv[])
{
  FILE *ff;
  char *month_name[100];
  int day_n,month,i,j;
  int month_day[12];
  float sig1[SLEN],sig2[SLEN];
  SAC_HD tysac1,tysac2;
  if(argc != 3)
    {
      printf("Usage: Calc_cor_time rec_file_1 rec_file_2\n");
      exit(-1);
    }
//  if((ff = fopen(argv[1], "r"))==NULL)
//  sprintf(fname1,"%s",argv[1]);
  if ( read_sac (argv[1], sig1, &tysac1, SLEN) == NULL ) {
    fprintf(stderr,"file %s not found\n", argv[1]);
    return 0;
  }
  if ( read_sac (argv[2], sig2, &tysac2, SLEN) == NULL ) {
    fprintf(stderr,"file %s not found\n", argv[2]);
    return 0;
  }
  i=0;
  month_day[0]=31;month_name[0]="JAN";
  month_day[1]=29;month_name[1]="FEB";
  month_day[2]=31;month_name[2]="MAR";
  month_day[3]=30;month_name[3]="APR";
  month_day[4]=31;month_name[4]="MAY";
  month_day[5]=30;month_name[5]="JUN";
  month_day[6]=31;month_name[6]="JUL";
  month_day[7]=31;month_name[7]="AUG";
  month_day[8]=30;month_name[8]="SEP";
  month_day[9]=31;month_name[9]="OCT";
  month_day[10]=30;month_name[10]="NOV";
  month_day[11]=31;month_name[11]="DEC";
   for(month=0;month<12;month++){
     day_n=0;
     for(j=1;j<=month_day[month];j++){
        if(sig1[i+j]*sig2[i+j]>0)day_n++;
//        printf("%d\n",day_n);
     }
     i+=month_day[month];
     printf("month: %s  num_of_days: %d\n",month_name[month],day_n);
   }
//  day_n=tysac.npts;
}
