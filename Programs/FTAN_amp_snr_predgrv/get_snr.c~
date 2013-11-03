#include <stdio.h>
#include <math.h>

#define SLEN 400000

void gaufilt_(double *alpha,double *c_per,
              double *dt,int *n, float seis_in[], float seis_out[]);

int get_snr(float *sei, int nsample, double dt, double dist, double b, double *c_per, double *g_vel, int nper, double *amp_max, double *snr2)
{
  double minT,maxT,signalmax,noisemax,noiserms;
  double alpha=20.,num,e;
  int k, i, ib, ie;
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
  ib = (int)floor(minT/dt);
  ie = (int)ceil(maxT/dt);

    signalmax=0;
    noisemax=0;
    for(i=ib;i<ie;i++) {
      if(seis_out[i] < 0) num = seis_out[i]*(-1);
      else num = seis_out[i];
      if(num>signalmax)
      signalmax=num;
    }
    noiserms=0;
    for(i=e-1000;i<e-500;i++)
      noiserms += pow(seis_out[i],2);
    noiserms=sqrt(noiserms/(500-1));
    snr2[k]=signalmax/noiserms;
    amp_max[k]=signalmax;
  }
  return 1;
}

