//This code extract the amplitude info from cross-correlation SAC file.
//Amplitude corrected with Stacking day number.
//Input file list with file name in the first colum. 
//Put day numbers in the second colum (optional) to check the day number info in the SAC files. 

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


float sig0[SLEN],seis_out[SLEN],sig1[SLEN],sig2[SLEN];
char fname[300];
SAC_HD tysac;
/*--------------------------------------------------------------*/
int get_snr(char *fname,int day_num_check,char *comp)
/*--------------------------------------------------------------*/
{
  FILE *f2;
  double dist,minV,maxV,minT,maxT,window,signalmax,noisemax,noiserms;
  double maxP=35.0;
  int nf=6;
  double per[nf],f[nf],f1[22],noise[nf],snr[nf],norm_amp_day[nf],norm_amp_noise[nf],amp_max[nf];
  double minP=5.0, num;
  double b,e,fb,fe,step;
//  double per1, per2, per3, per4;
  int k,i,ii,bi,ei,day_num;
//  char fname1[300], fltname[200], tstr[300];
  char *fname1, *fltname, *tstr;
  fname1=(char *) malloc(300*sizeof(char));
  fltname = (char *) malloc(300*sizeof(char));
  tstr = (char *) malloc(300*sizeof(char));

  printf("fname is %s\n", fname);

  /*---------------- reading sac file  -------------------*/
//system("pwd\n");
  if ( read_sac (fname, sig0, &tysac, SLEN) == NULL ) {
    fprintf(stderr,"file %s not found\n", fname);
    	
    free(fname1);
    free(fltname);
    free(tstr);


    return 0;
  }
  day_num=(int)tysac.user0;
  if( (day_num_check != day_num) && (day_num_check != -1) ) {
     fprintf(stderr,"Recorded stacking day number mismatck between input list and sac header for file %s.\n", fname);
     return 0;
    }
  if( day_num < 60 ) {
     fprintf(stderr,"Stacking time less than 60 days.\n");
     return 0;
    }

  fb=1.0/maxP;
  fe=1.0/minP;
  step=(log(fb)-log(fe))/(nf-1);
  for(k=0;k<nf;k++) {
    f[k]=exp(log(fe)+k*step);
    per[k]=1/f[k];
  }

  b=tysac.b;
  e=tysac.e;
//  e=-b;
  bi=(int)b;
  ei=(int)e;
//printf("b: %f  e: %f\n",b,e);
//printf("bi: %d  ei: %d\n",bi,ei);
for(i=-bi;i<=-bi+ei;i++) {
    sig2[i+bi]=(sig0[i]+sig0[-bi*2-i])/2;
    }
if(strcmp(comp,"sym")==0){
  for(i=-bi;i<=-bi+ei;i++) { 
    sig1[i+bi]=(sig0[i]+sig0[-bi*2-i])/2;
//if(sig0[i]*sig0[i]>3)printf("i: %d  i: %f  -i: %f  s: %f\n",i+bi,sig0[i],sig0[-bi*2-i],sig1[i+bi]);
    }
  }
else if(strcmp(comp,"pos")==0){
  for(i=0;i<=ei;i++) {
    sig1[i]=sig0[i-bi];
    }
  }
else {
  for(i=0;i<=-bi;i++) {
    sig1[i]=sig0[-bi-i];
    }
  }
  dist = tysac.dist;
//printf("%f\n",dist);  
  double dt;
  int nsamples,n,npow;
  dt = tysac.delta;
  nsamples = tysac.npts;
  npow=1;
  n=nsamples;

  minV = 2.0;
  maxV = 4.0;
  minT = dist/maxV-maxP/2;
  maxT = dist/minV+maxP;
  if(minT<0)
    minT=0;
  if(maxT>e)
    maxT=e;
  window=maxT-minT;
double ff1,ff2,ff3,ff4;
//for(k=0;k<=21;k++)
//  {
//     f1[k] = exp(log(0.4)-0.2*k);
//  }
//for(k=1;k<21;k++)
//  {
//      if(k>7 && k<17){continue;}
//      ff1 = f1[k+1]/1.25;
//      ff2 = f1[k+1];        
//      ff3 = f1[k-1];
//      ff4 = f1[k-1]*1.25;
//
//      filter4_(&ff1,&ff2,&ff3,&ff4,&npow,&dt,&n,sig1,seis_out);
  ff1=1./6./1.25;
  ff2=1./6.;
  ff3=1./4.;
  ff4=1./4.*1.25;
  filter4_(&ff1,&ff2,&ff3,&ff4,&npow,&dt,&n,sig2,seis_out);
  noiserms=0.;
  ii=0;
int noise_b=1500;
if(maxT+500>1500)noise_b=(int)maxT+500;
    for(i=noise_b;i<2500;i++) {
//      noise[k]=seis_out[i]*seis_out[i]+noise[k];
      noiserms=seis_out[i]*seis_out[i]+noiserms;
      ii++;
    }
    noiserms=sqrt(noiserms/(ii-1));
    //fprintf(stderr,"%lf\n",noiserms);
//  }
//  sprintf(fname1,"%s_amp_l_norm.txt",fname);
//  f2=fopen(fname1,"w");
//  for(i=1;i<21;i++) {
//    if(i>7 && i<17){continue;}
//    fprintf(f2,"%lf  normd_amp amp_max %lf\n",1/f1[i],noise[i]);
//  }
//  fclose(f2);
//  free(fname1);
//  return 1;

//  char buff[300];
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
      else { 

         free(fname1);
         free(fltname);
         free(tstr);
         return 1;

         }
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

    noise[k]=0.;
    noise_b=(int)maxT+500;
    for(i=noise_b;i<noise_b+1000;i++) {
      noise[k]=seis_out[i]*seis_out[i]+noise[k];
    }
    noise[k]=sqrt(noise[k]/(1000-1));

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

//printf("signalmax: %f  noiserms: %f\n",signalmax,noiserms);
    snr[k]=signalmax/noise[k];
    norm_amp_day[k]=signalmax/day_num;
//    norm_amp_noise[k]=signalmax/noiserms/noiserms;
//    norm_amp[k]=snr[k]/pow(noiserms,2./3.);
    amp_max[k]=signalmax;
//printf("period:%lf noiserms_square:%lf normd_amp:%lf amp:%lf\n",per[k],noiserms*noiserms,norm_amp[k],amp_max[k]);
    //    snr[k]=signalmax/noisemax;      
  }
//printf("%lf\n",noiserms*noiserms);
  sprintf(fname1,"%s_amp_l_norm_%s.txt",fname,comp);
  f2=fopen(fname1,"w");
  if(snr[(int)(nf/2)]<3 || snr[nf-1-(int)(nf/2)]<3){
    fprintf(stderr,"SNR of %s: %f. less than 3. skipped!\n",fname,snr[(int)(nf/2)] );
    sprintf(tstr,"rm -f %s\n",fname1);
    system(tstr);
    
    free(fname1);
    free(fltname);
    free(tstr);
    fclose(f2);
    return 0;
  }
  if(day_num<50){
//  if(noiserms<0.015){
    fprintf(stderr,"stacking time of %s: %d is too short. skipped!\n",fname,day_num );
    sprintf(tstr,"rm -f %s\n",fname1);
    system(tstr);

    free(fname1);
    free(fltname);
    free(tstr);
    fclose(f2);
    return 0;
  }
  for(i=1;i<nf-1;i++) {
    fprintf(f2,"%lf %g %lf %d %f\n",per[i],norm_amp_day[i],amp_max[i],day_num,snr[i]);
  }
//  printf("maxP=%f minP=%f maxV=%f minV=%f!!!\n",maxP,minP,maxV,minV);

    fclose(f2);

    free(fname1);
    free(fltname);
    free(tstr);

  return 1;
}


/*--------------------------------------------------------------*/
int main (int argc, char *argv[])
/*--------------------------------------------------------------*/
{
  
  FILE *ff;
  char filename[5000000];
  int day_num_check;
//  char filename1[300];
  if(argc != 3) 
    {
      printf("Usage: amp_length_normalize [list_file (optional: with_dayinfo)] [sym pos or neg]\n");
      exit(-1);
    }
  if((ff = fopen(argv[1], "r"))==NULL) {
    printf("Cannot open list.lst file. A file containing the names of");
    printf(" all cut files for getting spectral SNR is required to run spectral_snr*.\n");
    exit(1);
  }
  for(;;) {
     day_num_check=-1;
     if(fscanf(ff,"%s %d",&(filename[0]),&day_num_check)==EOF)
        break;
//    sprintf(filename1,"/utera/tianye/data_ASN_TA/%s",filename);
     get_snr(filename,day_num_check,argv[2]);
    }
fclose(ff);
return 1;
}
