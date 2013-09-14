/*
 * The sample of test driver for FTAN with phase match filter for
 * subroutines aftanpg and aftanipg
 */
#define MAIN
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

void gaufilt_(double *alpha,double *c_per,
              double *dt,int *n, float seis_in[], float seis_out[]);

void filter4_(double *f1,double *f2,double *f3,double *f4,
              double *dt,int *n, float seis_in[],
              int *ns,double *dom);

/*--------------------------------------------------------------*/
int get_snr(float *sei, int nsample, double dt, double dist, double b, double *c_per, double *g_vel, int nper, double *amp_max, double *snr2)
/*--------------------------------------------------------------*/
{
  FILE *fp1, *f2;
  double minT,maxT,window,signalmax,noisemax,noiserms;
//  int nf=64;
//  double f[nf];
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

/*--------------------------------------------------------------*/
int main (int argc, char *argv[])
{
  static int n, npoints, nfin, nfout1, nfout2, ierr, nprpv;
  static double t0, dt, delta, vmin, vmax, tmin, tmax;
  static double snr, tresh, ffact, perc, taperl,fmatch,piover4;
  static float sei[16384], sei_p[32768], sei_n[16384];
  static double arr1[100][8],arr2[100][7];
  static double c_per[100],g_vel[100],amp_p[100],amp_n[100];
  static double tamp, ampo[32][32768], pred[2][300];
  static int nrow, ncol, npred;
  static double prpvper[300],prpvvel[300]; // phase vel prediction files

  double snr_p[64], snr_n[64];
  double f1,f2,f3,f4,dom_am;
  char  *p,name[160],name1[160],buf[200],str[160],phvelname[160],root[160];
  char ev[20], sta2[6], amp_name[100];
  FILE  *in, *fd, *inv, *fas;
  int   i, j, flag, k, len, n_am, i_am;
  int   nn,sac = 1; // =1 - SAC, =0 - ftat files

// input command line arguments treatment
  if(argc != 3) {
      printf("Usage: aftan_amp [parameter file] [Model Path]\n");
      exit(-1);
  }
// open and read contents of parameter file
  if((in = fopen(argv[1],"r")) == NULL) {
      printf("Can not find file %s.\n",argv[1]);
      exit(1);
  }
  while((n = fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %d",
             &piover4,&vmin,&vmax,&tmin,&tmax,&tresh,&ffact,&taperl,&snr,&fmatch,
             name,&flag)) != EOF) { // start main loop
  strcpy(root,name);
  p = strrchr(root,'.');
  *(p+1) = '\0';
  //strcpy(phvelname,root);
  //strcat(phvelname,"SAC_PHP");

      if(n == 0 || n != 12) break;

      printf("vmin= %lf, vmax= %lf, tmin= %lf, tmax= %lf\n",
              vmin,vmax,tmin,tmax);
// remove quotes from file names
      j = 0;
      for(i = 0; i < strlen(name); i++) {
          if(name[i] == '\'' || name[i] == '\"') continue;
          name[j] = name[i]; j++;
      }
      name[j] = '\0';
      printf("Tresh= %lf, Filter factor= %lf, SNR= %lf, Match = %lf\nData file name=%s\n",
             tresh,ffact,snr,fmatch,name);
// if presents, read phase velocity prediction file
// ---
nprpv = -1;
                                                                                
/*
 *  read phase velocity information
 */
    strcpy(ev,strtok(root,"."));
    strcpy(sta2,strtok(NULL,"."));
//    sscanf(root, "%[A-Z,a-z,0-9]_%[A-Z,a-z,0-9].", sta1, sta2);

    sprintf(phvelname, "%s/%s/%s_%s.dat",argv[2],ev,ev,sta2);
   printf("%s\n",phvelname);
//    fprintf(stderr, "predicted phase velocity %s \n",phvelname);
  if((inv = fopen(phvelname,"r")) == NULL) {
      printf("Can not find model file %s. Skipped\n",phvelname);
      nprpv = 0;
      continue;
  }

  fgets(buf,200,inv);
  while(fgets(buf,200,inv) != NULL) {
         if(nprpv == -1) { nprpv++; continue; }
         if((n = sscanf(buf,"%lf %lf",&prpvper[nprpv],&prpvvel[nprpv])) < 2) break;
         nprpv++;
     }
         fclose(inv);
         printf("Phase velocity prediction file name= %s\n",phvelname);

//  next:
/*
 *   read SAC or ascii data 
 */
      readdata(sac,name,&n,&dt,&delta,&t0,sei_p);
      if(dt>1e10 || n<=1) continue;

      if(flag==1) {
         len=(n-1)/2;
         for(k=0;k<=len;k++) sei_n[k]=sei_p[len-k]; 
         for(k=0;k<=len;k++) sei_p[k]=sei_p[len+k];
         for(k=0;k<=len;k++) sei[k]=(sei_p[k]+sei_n[k])/2.;
         n=len+1;
	 t0 += len*dt;
        }
      else
         for(k=0;k<n;k++) sei[k]=sei_p[k];
/*
 * Read group velocity prediction file
 */
//      strcpy(name1,root);
//      strcat(name1,"SAC_GRP");
//      sscanf(root, "%[A-Z,a-z,0-9]_%[A-Z,a-z,0-9].", sta1, sta2);
      //sprintf(name1, "/home/linf/California/PRED_DISP/COR_%s_%s.SAC_PRED", sta1, sta2);
      //printf("Group velocity prediction curve: %s\n",name1);
      //if((fd = fopen(name1,"r")) == NULL) {
      //    printf("Can not find file %s.\n",name1);
      //    exit(1);
      //}
      //i = 0;
      //fgets(str,100,fd);
      //while((nn = fscanf(fd,"%lf %lf",&pred[0][i],&pred[1][i])) == 2) i++;
      //npred = i;
      //fclose(fd);
/* ====================================================================
 * Parameters for aftanipg function:
 * Input parameters:
 * piover4 - phase shift = pi/4*piover4, for cross-correlation
 *           piover4 should be -1.0 !!!!     (double)
 * n       - number of input samples, (int)
 * sei     - input array length of n, (float)
 * t0      - time shift of SAC file in seconds, (double)
 * dt      - sampling rate in seconds, (double)
 * delta   - distance, km (double)
 * vmin    - minimal group velocity, km/s (double)
 * vmax    - maximal value of the group velocity, km/s (double)
 * tmin    - minimal period, s (double)
 * tmax    - maximal period, s (double)
 * tresh   - treshold, usually = 10, (double)
 * ffact   - factor to automatic filter parameter, (double)
 * perc    - minimal length of of output segment vs freq. range, % (double)
 * npoints - max number points in jump, (int)
 * taperl  - factor for the left end seismogram tapering,
 *           taper = taperl*tmax,    (double)
 * nfin    - starting number of frequencies, nfin <= 32, (int)
 * snr     - phase match filter parameter, spectra ratio to
 *           determine cutting point    (double)
 * fmatch  - factor to length of phase matching window
 * npred   - length of prediction table
 * pred    - prediction table: pred[0][] - periods in sec,
 *                             pred[1][] - pedicted velocity, km/s
 * flag    - Input file type: 0 for single-sided, 1 for double-sided
 * ==========================================================
 * Output parameters are placed in 2-D arrays arr1 and arr2,
 * arr1 contains preliminary results and arr2 - final.
 * ==========================================================
 * nfout1 - output number of frequencies for arr1, (int)
 * arr1   - the first nfout1 raws contain preliminary data,
 *          (double arr1[n][5], n >= nfout1)
 *          arr1[:,0] -  central periods, s (double)
 *          arr1[:,1] -  apparent periods, s (double)
 *          arr1[:,2] -  group velocities, km/s (double)
 *          arr1[:,3] -  phase velocities, km/s (double)
 *          arr1[:,4] -  amplitudes, Db (double)
 *          arr1[:,5] -  discrimination function, (double)
 *          arr1[:,6] -  signal/noise ratio, Db (double)
 *          arr1[:,7] -  maximum half width, s (double)
 * nfout2 - output number of frequencies for arr2, (int)
 *          If nfout2 == 0, no final result.
 * arr2   - the first nfout2 raws contains final data,
 *          (double arr2[n][5], n >= nfout2)
 *          arr2[:,0] -  central periods, s (double)
 *          arr2[:,1] -  apparent periods, s (double)
 *          arr2[:,2] -  group velocities, km/s (double)
 *          arr2[:,3] -  amplitudes, Db (double)
 *          arr2[:,4] -  signal/noise ratio, Db (double)
 *          arr2[:,5] -  maximum half width, s (double)
 *          tamp      -  time to the beginning of ampo table, s (double)
 *          nrow      -  number of rows in array ampo, (int)
 *          ncol      -  number of columns in array ampo, (int)
 *          ampo      -  Ftan amplitude array, Db, (double [32][32768])
 * ierr   - completion status, =0 - O.K.,           (int)
 *                             =1 - some problems occures
 *                             =2 - no final results
 */

//  t0      = 0.0;
  nfin    = 32;
  npoints = 5;        // only 3 points in jump
  perc    = 50.0;     // 50 % for output segment
//  taperl  = 2.0;      // factor to the left end tapering
  printf("pi/4 = %5.1lf, t0 = %9.3lf\n",piover4,t0);
  printf("#filters= %d, Perc= %6.2f %s, npoints= %d, Taper factor= %6.2f\n",
          nfin,perc,"%",npoints,taperl);
/* Call aftanipg function, FTAN + prediction         */

  // printf("FTAN + prediction curve\n");

  //ffact =2.0;
  //aftanipg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh,
  //      &ffact,&perc,&npoints,&taperl,&nfin,&snr,&fmatch,&npred,pred,
  //      &nprpv,prpvper,prpvvel,
  //      &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);
  //printres(dt,nfout1,arr1,nfout2,arr2,tamp,nrow,ncol,ampo,ierr,name,"_P");
  //if(nfout2 == 0) continue;   // break aftan sequence
  //printf("Tamp = %9.3lf, nrow = %d, ncol = %d\n",tamp,nrow,ncol);

/* Pre-whiten and record the amp factor */
  f1=1./tmax/1.25;
  f2=1./tmax;
  f3=1./tmin;
  f4=1./tmin/1.25;
  filter4_(&f1,&f2,&f3,&f4,&dt,&n,sei,&n_am,&dom_am);
//  printf("%lf  %lf %lf\n",amp_rec[180],amp_rec[1000],amp_rec[1600]);
//  printf("%d  %lf\n",n_am/2+1,dom_am);
/* FTAN with phase match filter. First Iteration. */

  printf("FTAN - the first iteration\n");
  ffact =1.0;
//printf("pi: %lf  t0: %lf  dt: %lf  delta: %lf  vmin: %lf  vmax: %lf  tmin: %lf  tmax: %lf\n",piover4, t0, dt, delta, vmin, vmax, tmin, tmax);
  aftanpg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh,
        &ffact,&perc,&npoints,&taperl,&nfin,&snr,&nprpv,prpvper,prpvvel,
        &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);
  printres(dt,nfout1,arr1,nfout2,arr2,tamp,nrow,ncol,ampo,ierr,name,"_1");
  if( ierr == -1 ) {printf("Problematic input file. Skipped!\n"); continue; }
  if(nfout2 == 0 ) continue;   // break aftan sequence
  printf("Tamp = %9.3lf, nrow = %d, ncol = %d\n",tamp,nrow,ncol);

/* Make prediction based on the first iteration               */

  npred = nfout2;
  tmin = arr2[0][1];
  tmax = arr2[nfout2-1][1];
  for(i = 0; i < nfout2; i++) {
      pred[0][i] = arr2[i][1];   // apparent period // central periods
      pred[1][i] = arr2[i][2];   // group velocities
  }

/* FTAN with phase with phase match filter. Second Iteration. */

  printf("FTAN - the second iteration (phase match filter)\n");
  ffact = 2.0;
  aftanipg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh,
        &ffact,&perc,&npoints,&taperl,&nfin,&snr,&fmatch,&npred,pred,
        &nprpv,prpvper,prpvvel,
        &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);
  printf("Tamp = %9.3lf, nrow = %d, ncol = %d\n",tamp,nrow,ncol);
  printres(dt,nfout1,arr1,nfout2,arr2,tamp,nrow,ncol,ampo,ierr,name,"_2");
  for(i = 0; i < nfout2; i++) {
      c_per[i]=arr2[i][0];
      g_vel[i]=arr2[i][2];
     }
  get_snr(sei_p,n,dt,delta,t0,c_per,g_vel,nfout2,amp_p,snr_p);
  if(flag==1) get_snr(sei_n,n,dt,delta,t0,c_per,g_vel,nfout2,amp_n,snr_n);
  sprintf(amp_name,"%s_amp_snr",name);
  if((fas=fopen(amp_name,"w"))==NULL) {
     printf("Cannot open file %s to write!\n");
     exit (1);
    }
  if(flag==1) for(i = 0; i < nfout2; i++)
     fprintf(fas,"%8.4f   %.5g  %8.4f  %.5g  %8.4f\n",arr2[i][1],amp_p[i],snr_p[i],amp_n[i],snr_n[i]);
  else for(i = 0; i < nfout2; i++) 
     fprintf(fas,"%8.4f   %.5g  %8.4f\n",arr2[i][1],amp_p[i],snr_p[i]);
  fclose(fas);
  }
  fclose(in);
  return 0;
}
