#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void filter4_(double *f1, double *f2, double *f3, double *f4, double *dt, int *n, float seis_in[], int *ns,double *dom);

int get_P(float *sei, int nsample, double dt, double dist, double b, double minT, double maxT, double *tvt, double *amp, double *snr)
/*--------------------------------------------------------------*/
{
  FILE *fp1, *f2;
  double window,signalmax,noisemax,noiserms;
  double alpha=20.,num,e;
  double dom;
  int i, ns, wb, we;
  char fname1[300], fltname[200];
  float seis_out[SLEN];

  wb=(int)(minT/dt+0.5);
  we=(int)(maxT/dt+0.5);
  e=b+(nsample-1)*dt;

  double f2=0.005, f1=f2/1.25, f3=0.3, f4=f3*1.25;
  filter4_(&f1,&f2,&f3,&f4,&dt,&n,seis_in,&ns,&dom);

  signalmax=0;
  noisemax=0;
  for(i=wb;i<we;i++) {
     if(seis_out[i] < 0) num = seis_out[i]*(-1);
     else num = seis_out[i];
     if(num>signalmax)
     signalmax=num;
  }
  noiserms=0;
  for(i=e-1000;i<e-500;i++) {
     noiserms=seis_out[i]*seis_out[i]+noiserms;
  }
  noiserms=sqrt(noiserms/(500-1));
  snr2[k]=signalmax/noiserms;
  amp[k]=signalmax;

  return 1;
}
