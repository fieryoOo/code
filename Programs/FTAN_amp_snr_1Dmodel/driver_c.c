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
#include "mysac64.h"
#include "/home/tianye/code/Programs/head/koftan.h"
#include "/home/tianye/code/Programs/head/gl_const.h"
#include "/home/tianye/code/Programs/head/mymacro.h"

#define SLEN 400000

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD);

void write_sac (char *fname, float *sig, SAC_HD *SHD);

/*--------------------------------------------------------------*/
int get_snr(float *sei, int nsample, double dt, double dist, double b, double *c_per, double *g_vel, int nper, double *amp_max, double *snr2)
/*--------------------------------------------------------------*/
{
  FILE *fp1, *f2;
  double minT,maxT,window,signalmax,noisemax,noiserms;
  int nf=64;
  double f[nf];
  double alpha=20.,num,e;
  int k,i, j,ii, ib, ie;
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

/*--------------------------------------------------------------*/
int pflag;
int main (int argc, char *argv[])
{
  static int n, npoints, nfin, nfout1, nfout2, ierr, nprpv;
  static double t0, dt, delta, vmin, vmax, tmin, tmax;
  static double snr, tresh, ffact, perc, taperl,fmatch,piover4;
  static float sei[16384], seiout[16384], sei_n[16384], *sei_p;
  static double arr1[100][8],arr2[100][7];
  static double c_per[100],g_vel[100],amp_p[100],amp_n[100];
  static double tamp, ampo[32][32768];
  static int nrow, ncol;
  static double prpvper[300],prpvvel[300]; // phase vel prediction files

  double snr_p[64], snr_n[64];
  double f1,f2,f3,f4,dom_am;
  char  *p,name[160],name1[160],buf[200],str[160],phvelname[160],root[160];
  char sta1[6], sta2[6], amp_name[100];
  FILE  *in, *fd, *inv, *fas;
  int   i, j, flag, k, len, n_am, i_am;
  int   nn,sac = 1; // =1 - SAC, =0 - ftat files
  SAC_HD shd;

// input command line arguments treatment
  if(argc!=3 && argc!=4 ) {
      printf("Usage: aftan_amp [parameter file] [file_pred_phvel] [out_flag(optional)]\n");
      exit(-1);
  }
/*---------------- out_flag --------------------
controls what files to output
0(default):     all the files
1:		only _1_DISP.1 and _amp_snr
2:              only _2_DISP.1 and _amp_snr
notice that amp_snrs are always measured based on
 the final dispersions (_2_DISP.1)
-----------------------------------------------*/
   pflag = 0;
   if(argc==4) pflag = atof(argv[3]);
   if(pflag<0 || pflag>2) {
      printf("Unknow out_flag: %d\n", pflag);
      exit(-1);
   }
// open and read contents of parameter file
  if((in = fopen(argv[1],"r")) == NULL) {
      printf("Can not find file %s.\n",argv[1]);
      exit(1);
  }
// read in 1D model
  if((inv = fopen(argv[2],"r")) == NULL) {
     printf("Cannot open model file %s.\n", argv[2]);
     exit(0);
  }
  nprpv = 0;
  while(fgets(buf,200,inv) != NULL) {
         if((n = sscanf(buf,"%lf %lf",&prpvper[nprpv],&prpvvel[nprpv])) < 2) break;
         nprpv++;
     }
  fclose(inv);
//  printf("Phase velocity prediction file name= %s\n",phvelname);
// main loop; read in one sac record at a time and do FTAN
  while((n = fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %d",
             &piover4,&vmin,&vmax,&tmin,&tmax,&tresh,&ffact,&taperl,&snr,&fmatch,
             name,&flag)) != EOF) {

      if(n == 0 || n != 12) break;

      //printf("vmin= %lf, vmax= %lf, tmin= %lf, tmax= %lf\n",
      //        vmin,vmax,tmin,tmax);
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
                                                                                
/*
 *  read phase velocity information
 */
//    strtok(root,"_");
//    strcpy(sta1,strtok(NULL,"_"));
//    strcpy(sta2,strtok(NULL,"."));
//    sscanf(root, "%[A-Z,a-z,0-9]_%[A-Z,a-z,0-9].", sta1, sta2);

//    sprintf(phvelname, "%s/%s/%s_%s.dat",argv[2],sta1,sta1,sta2);
//    fprintf(stderr, "predicted phase velocity %s \n",phvelname);
/*
  if((inv = fopen(phvelname,"r")) == NULL) {
      printf("Can not find file %s. Use the inversed path instead\n",phvelname);
      sprintf(phvelname, "%s/%s/%s_%s.dat",argv[2],sta2,sta2,sta1);
      if((inv = fopen(phvelname,"r")) == NULL) {
         printf("Can not find file %s either.\n",phvelname);
         nprpv = 0;
         //goto next;
         continue;
        }
     }
*/
//  next:
/*
 *   read SAC or ascii data 
 */
//      readdata(sac,name,&n,&dt,&delta,&t0,sei_p);
   if( read_sac(name, &sei_p, &shd) == NULL ) {
      fprintf(stderr, "ERROR(read_sac): %s\n", name);
      exit(-1);
   }
   n = shd.npts; dt = shd.delta;
   delta = shd.dist; t0 = shd.b;

      if(flag==1) { //the input sac is marked as double-sided. compute symetric compenent
         len=(n-1)/2;
         for(k=0;k<=len;k++) sei_n[k]=sei_p[len-k]; 
         for(k=0;k<=len;k++) sei_p[k]=sei_p[len+k];
         for(k=0;k<=len;k++) sei[k]=(sei_p[k]+sei_n[k])/2.;
         n=len+1;
	 t0 += len*dt;
      }
      else for(k=0;k<n;k++) sei[k]=sei_p[k]; //single-sided. copy
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
  while((n = fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %d",
             &piover4,&vmin,&vmax,&tmin,&tmax,&tresh,&ffact,&taperl,&snr,&fmatch,
             name,&flag)) != EOF) {
----------------------------------------------------------------------
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
 *          ampo      -  Ftan amplitude array, Db, (double [32][32768])
 * ierr   - completion status, =0 - O.K.,           (int)
 *                             =1 - some problems occures
 *                             =2 - no final results
 */

//  t0      = 0.0;
  nfin    = 32;
  npoints = 10;        // only 3 points in jump
  perc    = 35.0;     // output if the percentage of measurable frequecy range is greater than 35%
//  taperl  = 2.0;      // factor to the left end tapering
  //printf("pi/4 = %5.1lf, t0 = %9.3lf\n",piover4,t0);
  //printf("#filters= %d, Perc= %6.2f %s, npoints= %d, Taper factor= %6.2f\n",
  //        nfin,perc,"%",npoints,taperl);
/* Call aftanipg function, FTAN + prediction         */

  // printf("FTAN + prediction curve\n");

  //ffact =2.0;
  //aftanipg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh,
  //      &ffact,&perc,&npoints,&taperl,&nfin,&snr,&fmatch,&npred,pred,
  //      &nprpv,prpvper,prpvvel,
  //      &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);
  //printres(dt,nfout1,arr1,nfout2,arr2,tamp,nrow,ncol,ampo,ierr,name,"_P",delta);
  //if(nfout2 == 0) continue;   // break aftan sequence
  //printf("Tamp = %9.3lf, nrow = %d, ncol = %d\n",tamp,nrow,ncol);

/* Pre-whiten and record the amp factor */
  f1=1./tmax/1.25;
  f2=1./tmax;
  f3=1./tmin;
  f4=1./tmin/1.25;
//  filter4_(&f1,&f2,&f3,&f4,&dt,&n,sei,&n_am,&dom_am);
//  printf("%lf  %lf %lf\n",amp_rec[180],amp_rec[1000],amp_rec[1600]);
//  printf("%d  %lf\n",n_am/2+1,dom_am);
/* FTAN (without?) phase match filter. First Iteration. */

  printf("FTAN - the first iteration\n");
  double tresh1 = tresh * 1.;
  aftanpg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh1,
        &ffact,&perc,&npoints,&taperl,&nfin,&snr,&nprpv,prpvper,prpvvel,
        &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);
  if( pflag==0 || pflag==1 ) printres(dt,nfout1,arr1,nfout2,arr2,tamp,nrow,ncol,ampo,ierr,name,"_1",delta);
  if(nfout2 == 0) continue;   // break aftan sequence
  printf("Tamp = %9.3lf, nrow = %d, ncol = %d\n",tamp,nrow,ncol);

/* Read in the predicted group dispersion. (or make prediction based on the first iteration.) */
  static double pred[2][300];
  static int npred;

  npred = nfout2;
  tmin = arr2[0][1];
  tmax = arr2[nfout2-1][1];
  for(i = 0; i < nfout2; i++) {
      pred[0][i] = arr2[i][1];   // apparent periods
      pred[1][i] = arr2[i][2];   // group velocities
  }

/* FTAN with phase match filter. Second Iteration. */

  printf("FTAN - the second iteration (phase match filter)\n");
  int cuttype = 0;
  aftanipg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh,
        &ffact,&perc,&npoints,&taperl,&nfin,&snr,&fmatch,&npred,pred,
        &cuttype,&nprpv,prpvper,prpvvel,seiout,
        &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);
  printf("Tamp = %9.3lf, nrow = %d, ncol = %d\n",tamp,nrow,ncol);
  if( pflag==0 || pflag==2 ) printres(dt,nfout1,arr1,nfout2,arr2,tamp,nrow,ncol,ampo,ierr,name,"_2",delta);
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
