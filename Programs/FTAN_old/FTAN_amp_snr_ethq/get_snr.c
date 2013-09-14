#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aftan.h"
#include "/home/tianye/code/Programs/head/mysac64.h"
#include "/home/tianye/code/Programs/head/sac_db64.h"
#include "/home/tianye/code/Programs/head/koftan.h"
#include "/home/tianye/code/Programs/head/gl_const.h"
#include "/home/tianye/code/Programs/head/mymacro.h"
#define SLEN 400000

int get_snr(float *sei, int nsample, double dt, double dist, double b, double *c_per, double *g_vel, int nper, double *amp_max, double *snr2)
{
  FILE *fp1, *f2;
  double minT,maxT,window,signalmax,noisemax,noiserms;
  double alpha=20.,num,e;
  int k,i, j,ii, iwin, ib, ie;
  char fname1[300], fltname[200];
  float seis_out[SLEN];

  e=b+(nsample-1)*dt;
//  printf ("nsig: %d  dt: %f  dist: %f  b: %f  e: %f  nper: %d  alpha: %f\n",nsample,dt,dist,b,e,nper,alpha);
  for(k = 0; k < nper; k++) {
//     printf("centre_per: %8.4lf   Group_vel: %6.4lf\n",c_per[k],g_vel[k]);

  gaufilt_(&alpha, &c_per[k],&dt,&nsample,sei,seis_out);
  minT = dist/g_vel[k]-c_per[k]/2.;
  maxT = dist/g_vel[k]+c_per[k]/2.;
  if(minT<b)
    minT=b;
  if(maxT>e)
    maxT=e;
  window=maxT-minT;

    signalmax=0;
    noisemax=0;
    for(i=(int)minT;i<maxT;i++) {
      if(seis_out[i] < 0) num = seis_out[i]*(-1);
      else num = seis_out[i];
      if(num>signalmax) {
         signalmax=num;
         ib=i;
      }
    }
    if(e-ib<700) {
       printf("time series not long enough for computing snr!\n");
       snr2[k]=0;
       amp_max[k]=signalmax;
       continue;
    }
    else if(e-ib>=700 && e-ib<1100) {
       ie=e-100; ib+=500; iwin=ie-ib; }
    else { ib+=500; iwin=499; ie=ib+iwin; }
    noiserms=0;
    for(i=ib;i<ie;i++) {
      noiserms=seis_out[i]*seis_out[i]+noiserms;
    }
    noiserms=sqrt(noiserms/iwin);
    snr2[k]=signalmax/noiserms;
    amp_max[k]=signalmax;
//    printf("amp: %g  snr: %f\n",amp_max[k],snr2[k]);
  }
  return 1;
}

