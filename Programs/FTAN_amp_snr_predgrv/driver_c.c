/*
 * The sample of test driver for FTAN with phase match filter for
 * subroutines aftanpg and aftanipg
 */
/* ====================================================================
 * Parameters for aftanipg function:
----------------------------------------------------------------------
  while((n = fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %d",
             &piover4,&vmin,&vmax,&tmin,&tmax,&tresh,&ffact,&taperl,&snr,&fmatch,
             name,&flag)) != EOF) {
----------------------------------------------------------------------
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
 * perc    - minimal length of output segment vs freq. range, % (double)
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
 *          ampo      -  Ftan amplitude array, Db, (double [32][2*NMAX])
 * ierr   - completion status, =0 - O.K.,           (int)
 *                             =1 - some problems occures
 *                             =2 - no final results
 */


#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aftan.h"
#include "mysac64.h"
#include "/home/tianye/code/Programs/head/koftan.h"
#include "/home/tianye/code/Programs/head/gl_const.h"
#include "/home/tianye/code/Programs/head/mymacro.h"

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD);

void write_sac (char *fname, float *sig, SAC_HD *SHD);

int get_snr(float *sei, int nsample, double dt, double dist, double b, double *c_per, double *g_vel, int nper, double *amp_max, double *snr2);


/*--------------------------------------------------------------*/
#define NMAX 16384
int pflag = 1;
int main (int argc, char *argv[])
{
  static int n, npoints, nfin, nfout1, nfout2, ierr, nprpv;
  static double t0, dt, delta, vmin, vmax, tmin, tmax;
  static double snr, tresh, ffact, perc, taperl,fmatch,piover4;
  static double snr2, tresh2, ffact2, taperl2,fmatch2;
  static float sei[NMAX], seiout[NMAX], sei_n[NMAX], *sei_p;
  static double arr1[100][8],arr2[100][7];
  static double c_per[100],g_vel[100],amp_p[100],amp_n[100];
  static double tamp, ampo[32][NMAX*2];
  static int nrow, ncol;
  static double prpvper[300],prpvvel[300]; // phase vel prediction files

  double snr_p[64], snr_n[64];
  double f1,f2,f3,f4,dom_am;
  char name[150],buff[300];
  char amp_name[160];
  FILE *in, *inv, *fas;
  int i, j, flag, k, len, n_am;
  SAC_HD shd;

// input command line arguments treatment
  if(argc!=4 && argc!=5 ) {
      printf("Usage: aftan_amp [parameter file] [file_pred_phvel] [file_pred_grvel] [out_flag(optional)]\n");
      exit(-1);
  }
/*---------------- out_flag --------------------
controls what files to output
1(default):	all the files
0:		only _2_DISP.1 and _amp_snr
-----------------------------------------------*/
   pflag = 1;
   if(argc==5) pflag = atof(argv[4]);
   if(pflag!=0 && pflag!=1) {
      printf("Unknow out_flag: %d\n", pflag);
      exit(-1);
   }
// read in 1D model
  if((inv = fopen(argv[2],"r")) == NULL) {
     printf("Cannot open model file %s.\n", argv[2]);
     exit(0);
  }
  nprpv = 0;
  while(fgets(buff,300,inv) != NULL) {
         if((n = sscanf(buff,"%lf %lf",&prpvper[nprpv],&prpvvel[nprpv])) < 2) break;
         nprpv++;
     }
  fclose(inv);
// main loop
  // open and read contents of parameter file
  if((in = fopen(argv[1],"r")) == NULL) {
      printf("Can not find file %s.\n",argv[1]);
      exit(1);
  }
  while((n = fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %d",
             &piover4,&vmin,&vmax,&tmin,&tmax,&tresh,&ffact,&taperl,&snr,&fmatch, &tresh2,
	     &ffact2,&taperl2,&snr2,&fmatch2,name,&flag)) != EOF) {

      if(n == 0 || n != 17) break;

  // remove quotes from file names
      j = 0;
      for(i = 0; i < strlen(name); i++) {
          if(name[i] == '\'' || name[i] == '\"') continue;
          name[j] = name[i]; j++;
      }
      name[j] = '\0';

      printf("Tresh= %lf, Filter factor= %lf, SNR= %lf, Match = %lf\nData file name=%s\n",
             tresh,ffact,snr,fmatch,name);
  // read SAC or ascii data 
   int sac = 1; // =1 - SAC, =0 - ftat files
   //sei_p = (float *) malloc (2*NMAX * sizeof(float));
   //readdata(sac,name,&n,&dt,&delta,&t0,sei_p);
   if( read_sac(name, &sei_p, &shd) == NULL ) {
      fprintf(stderr, "ERROR(read_sac): %s\n", name);
      exit(-1);
   }
   n = shd.npts; dt = shd.delta;
   delta = shd.dist; t0 = shd.b;
   
  // prepare single sided record for FTAN
      if(flag==1) { //the input sac is marked as double-sided. compute symetric compenent
	 len = shd.npts/2;
	 if( len > NMAX ) {
	    fprintf(stderr, "ERROR(NMAX): len=%d > NMAX=%d!\n", len, NMAX);
	    exit(-1);
	 }
         for(k=0;k<=len;k++) sei_n[k]=sei_p[len-k]; 
         for(k=0;k<=len;k++) sei_p[k]=sei_p[len+k];
         for(k=0;k<=len;k++) sei[k]=(sei_p[k]+sei_n[k])/2.;
         n = len+1; t0 += len*dt;
	 shd.npts = n; shd.b = t0;
      }
      else {
	 if( n > NMAX ) {
	    fprintf(stderr, "ERROR(NMAX): n=%d > NMAX=%d!\n", n, NMAX);
	    exit(-1);
	 }
	 for(k=0;k<n;k++) sei[k]=sei_p[k]; //single-sided. copy
      }

  // Pre-whiten and record the amp factor
  f1=1./tmax/1.25;
  f2=1./tmax;
  f3=1./tmin;
  f4=1./tmin/1.25;
  //filter4_(&f1,&f2,&f3,&f4,&dt,&n,sei,&n_am,&dom_am);

  // set parameters
  nfin    = 32;
  npoints = 10;        // only 3 points in jump
//  taperl  = 2.0;      // factor to the left end tapering

/* FTAN without phase match filter to construct the original FTAN diagram. First Iteration. */
  perc    = 40.0; // output if the percentage of measurable frequecy range is greater than 10%
  printf("FTAN - the first iteration\n");
  double tresh1 = tresh * 1.;
  aftanpg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh1,
        &ffact,&perc,&npoints,&taperl,&nfin,&snr,&nprpv,prpvper,prpvvel,
        &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);
  if(pflag) printres(dt,nfout1,arr1,nfout2,arr2,tamp,nrow,ncol,ampo,ierr,name,"_1",delta);
  if(nfout2 == 0) continue;   // break aftan sequence 
  //printf("Tamp = %9.3lf, nrow = %d, ncol = %d\n",tamp,nrow,ncol);

/* Read in the predicted group dispersion. (or make prediction based on the first iteration.) */
  static double pred[2][300];
  static int npred;
  if( (inv=fopen(argv[3], "r")) == NULL ) {
     perror("Failed to open file ");
     exit(-1);
  }
  npred = 0; 
  float tminp = 9999., tmaxp = -1.;
  while( fgets(buff,300,inv) != NULL ) {
     if( sscanf(buff,"%lf %lf",&pred[0][npred],&pred[1][npred]) < 2) break;
     if(pred[1][npred] != pred[1][npred]) continue;
     if( pred[0][npred] < tmin || pred[0][npred] > tmax ) continue;
     if(tminp>pred[0][npred]) tminp = pred[0][npred];
     if(tmaxp<pred[0][npred]) tmaxp = pred[0][npred];
    // printf("%d %lf %lf\n", npred, pred[0][npred], pred[1][npred]);
     npred++;
  }
  fclose(inv);
  //printf("%f %f  %f %f\n", tmin, tmax, tminp, tmaxp); exit(0);
  //if( tminp > tmin ) tmin = tminp;
  //if( tmaxp < tmax ) tmax = tmaxp;
  tmin = tminp; tmax = tmaxp;

/* Pre-whiten and record the amp factor */
  f1=1./(tmax*1.25);
  f2=1./tmax;
  f3=1./tmin;
  f4=1./tmin*1.25;
  float amprec[NMAX], ampavg;
  filter4_(&f1,&f2,&f3,&f4,&dt,&n,sei,&n_am,&dom_am,amprec);
/*
  int ib = (int)ceil(f1/dom_am);
  for(i=ib;i<n_am;i++) {
     if(dom_am*i >= f4) break;
     ampavg += amprec[i];
  }
  ampavg = (i-ib)/ampavg;
*/
/* FTAN with (1st, wide) phase match filter to trach energy around prediction. Second Iteration.
 cuttype = 1 tells tgauss() to cut the anti-dispersed diagram by 
 strictly following the group prediction instead of at the max energy */
  int cuttype = 1;
  perc = 40.0;
  printf("FTAN - the second iteration (phase match filter)\n");
  aftanipg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh, // 11 params
        &ffact,&perc,&npoints,&taperl,&nfin,&snr,&fmatch,&npred,pred,	   // 9
        &cuttype,&nprpv,prpvper,prpvvel,seiout,				   // 4
        &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);	   // 9
  //printres(dt,nfout1,arr1,nfout2,arr2,tamp,nrow,ncol,ampo,ierr,name,"_2",delta);
  //sprintf(amp_name, "%s_cld", name);
  //write_sac (amp_name, seiout, &shd);

/* FTAN with (2nd, narrow) phase match filter. Third Iteration */
  npred = nfout2;
  tmin = arr2[0][1];
  tmax = arr2[nfout2-1][1];
  pred[0][0] = arr2[0][1]; pred[1][0] = arr2[0][2];
  for(i=1,j=1; i<nfout2; i++) {
      pred[0][j] = arr2[i][1];   // apparent periods
      if( pred[0][j] <= pred[0][j-1] ) continue;
      pred[1][j] = arr2[i][2];   // group velocities
      //printf("%d %f %f\n", j, pred[0][j], pred[1][j]); //aa
      j++;
  }
  cuttype = 0;
  perc = 40.0;
  printf("FTAN - the third iteration (phase match filter)\n");
  printf("%f %f\n", tmin, tmax);
  aftanipg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh2, // 11 params
        &ffact2,&perc,&npoints,&taperl2,&nfin,&snr2,&fmatch2,&npred,pred,      // 9
        &cuttype,&nprpv,prpvper,prpvvel,seiout,                                   // 4
        &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);           // 9
  printres(dt,nfout1,arr1,nfout2,arr2,tamp,nrow,ncol,ampo,ierr,name,"_2",delta);
  sprintf(amp_name, "%s_cld", name);
  write_sac (amp_name, seiout, &shd);

/* comput amplitude and SNR based on the phase-match filter results */
  for(i = 0; i < nfout2; i++) {
      c_per[i]=arr2[i][0];
      g_vel[i]=arr2[i][2];
     }
  get_snr(sei_p,n,dt,delta,t0,c_per,g_vel,nfout2,amp_p,snr_p); // positive lag
  if(flag==1) get_snr(sei_n,n,dt,delta,t0,c_per,g_vel,nfout2,amp_n,snr_n); //negative lag if exists
  sprintf(amp_name,"%s_amp_snr",name);
  if((fas=fopen(amp_name,"w"))==NULL) {
     printf("Cannot open file %s to write!\n", amp_name);
     exit (1);
    }
  if(flag==1) for(i = 0; i < nfout2; i++)
     fprintf(fas,"%8.4f   %.5g  %8.4f  %.5g  %8.4f\n",arr2[i][1],amp_p[i],snr_p[i],amp_n[i],snr_n[i]);
  else for(i = 0; i < nfout2; i++) 
     fprintf(fas,"%8.4f   %.5g  %8.4f\n",arr2[i][1],amp_p[i],snr_p[i]);
  fclose(fas);
  }
  fclose(in);
  free(sei_p);
  return 0;
}
