#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "/home/tianye/code/Programs/head/mysac64.h"
#include "/home/tianye/code/Programs/head/sac_db64.h"

#define SLEN 400000
#define NLOC 20

int readdata(int sac,char *name,int *n,double *dt,double *delta, double *t0,float *sei);

void gaufilt_(double *alpha, double *c_per, double *dt, int *n, float seis_in[], float seis_out[]);

int get_P(float *sei, int nsample, double dt, double dist, double b, float wb, float we, double *tvt, double *amp, double *snr);

int main (int argc, char *argv[])
{
  FILE *ff;
  static double t0, dt, delta,
  static float sei_p[32768];
  char Z_file[100], R_file[100], buff[300];
  int i, j, sac=1, n;
  double tvtZ[NLOC], ampZ[NLOC], snrZ[NLOC], tvtR[NLOC], ampR[NLOC], snrR[NLOC];
  double iinc;
  
  if((ff=fopen())==NULL) {
     cout<<""<<endl;
     return 0;
  }
  for(;;) {
     if((fgets(buff, 300, ff))==NULL) break;
     sscanf(buff, "%s %s %lf", &Z_file, &R_file, &iinc);

     readdata(sac,Z_file,&n,&dt,&delta,&t0,sei_p);
     get_P(sei_p, n, dt, delta, t0, wb, we, &tvtZ, &ampZ, &snrZ);

     readdata(sac,R_file,&n,&dt,&delta,&t0,sei_p);
     get_P(sei_p, n, dt, delta, t0, wb, we, &tvtR, &ampR, &snrR);

     iinc
     for(i=0;i<NLOC;i++)
        for(j=0;j<NLOC;j++) {
           if(fabs(tvtZ[i]-tvtE[j])>2) continue;
        }

  }
  fclose(ff);

  return 1;

}
