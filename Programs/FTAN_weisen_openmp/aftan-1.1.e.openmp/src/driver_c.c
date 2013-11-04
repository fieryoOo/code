/*
 * The sample of test driver for FTAN with phase match filter for
 * subroutines aftanpg and aftanipg
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aftan.h"
#include <omp.h>

//#include <iostream>

//using namespace std;
/*
#include "/home/weisen/PROGS_64/CC/head/mysac64.h"
#include "/home/weisen/PROGS_64/CC/head/sac_db64.h"
#include "/home/weisen/PROGS_64/CC/head/mymacro.h"
#include "/home/weisen/PROGS_64/CC/head/gl_const.h"
#include "/home/weisen/PROGS_64/CC/head/koftan.h"
*/

///home/wshen/my_code/CC/head

#include "/home/wshen/my_code/CC/head/mysac64.h"
#include "/home/wshen/my_code/CC/head/sac_db64.h"
#include "/home/wshen/my_code/CC/head/mymacro.h"
#include "/home/wshen/my_code/CC/head/gl_const.h"
#include "/home/wshen/my_code/CC/head/koftan.h"
//using namespace std;

#define SLEN 400000



/*--------------------------------------------------------------*/
int get_snr(float *sei, int nsample, double dt, double dist, double b, double *c_per, double *g_vel, int nper, double *amp_max, double *snr2, double *snrpc)
/*--------------------------------------------------------------*/
{
  FILE *fp1, *f2;
  double minT,maxT,window,signalmax,noisemax,noiserms,noiseprec,cvt;
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
  //fprintf(stderr,"cperk %g %d\n",c_per[k],nsample);
  gaufilt_(&alpha, &c_per[k],&dt,&nsample,sei,seis_out);
  minT = dist/g_vel[k]-c_per[k]/2.;
  maxT = dist/g_vel[k]+c_per[k]/2.;
  if(minT<b)
    minT=b;
  if(maxT>e)
    maxT=e;
  ib = (int)floor(minT/dt);
  ie = (int)ceil(maxT/dt);
  //fprintf(stderr,"ib: %g ie: %g dist: %g g_vel: %g\n",ib,ie,dist,g_vel[k]);

    signalmax=0;
    noisemax=0;
    for(i=ib;i<ie;i++) {
      if(seis_out[i] < 0) num = seis_out[i]*(-1);
      else num = seis_out[i];
      //fprintf(stderr,"num: %g\m",num);
      if(num>signalmax)
      signalmax=num;
    }
    noiserms=0;
    for(i=e-1000;i<e-500;i++)
      noiserms += pow(seis_out[i],2);
    noiserms=sqrt(noiserms/(500-1));
    snr2[k]=signalmax/noiserms;
    ii = 0;

    noiseprec = 0.;
    if (1000. < minT-6*c_per[k]) cvt = 1000.;
    else cvt = minT-6*c_per[k];

    for (i=0;i<cvt;i++) {
      noiseprec += pow(seis_out[i],2);
      ii ++;
      }
    if (ii> 300) {
      noiseprec = sqrt(noiseprec/(ii-1));
      snrpc[k] = signalmax/noiseprec;
      }
    else snrpc[k] = signalmax*-1.;
    amp_max[k]=signalmax;
  }
  return 1;
}



int main (int argc, char *argv[])
{
 /*
 static  int n, npoints, nfin, nfout1, nfout2, ierr, nprpv;
 static  double t0, dt, delta, vmin, vmax, tmin, tmax;
 static  double snr, tresh, ffact, perc, taperl,fmatch,piover4;
 static  float sei[32768],sei_in[32768];
 static  double arr1[100][8],arr2[100][7];
 static  double c_per[100],g_vel[100],amp_p[100],snr_p[100],snr_pc[100];
 static  double tamp, ampo[100][32768], pred[2][300];
 static  int nrow, ncol, npred;
 static  double prpvper[300],prpvvel[300]; */ /* phase vel prediction files  */

  /*
  char  *p,name[160],buf[200],phvelname[300],root[160],amp_name[160];
  FILE  *in, *inv, *fas;
  int   i,ii,pflag;
  int   sac = 1;  */ /* =1 - SAC, =0 - ftat files    */
  FILE *in;
  int j,n,ii;
/* input command line arguments treatment   */
  if(argc != 2) {
      printf("Usage: aftan4_c_test parameter_file\n");
      exit(-1);
  }
/* open and read contents of parameter file  */
  if((in = fopen(argv[1],"r")) == NULL) {
      printf("Can not find file %s\n",argv[1]);
      exit(1);
  }

  int nin;
  int NNN = 3000;
  double piover4s[NNN],vmins[NNN],vmaxs[NNN],tmins[NNN],tmaxs[NNN];
  double treshs[NNN],ffacts[NNN],taperls[NNN],snrs[NNN],fmatchs[NNN];
  int pflags[NNN];
  char names[NNN][160],phvelnames[NNN][300];
  
  j=0;
  while((n = fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %s %s",
            &(piover4s[j]),&vmins[j],&vmaxs[j],&tmins[j],&tmaxs[j],&treshs[j],&ffacts[j],
            &taperls[j],&snrs[j],&fmatchs[j],&pflags[j],
            names[j],phvelnames[j])) != EOF) { /* start main loop      */
  /*
  strcpy(root,name);
  p = strrchr(root,'.');
  *(p+1) = '\0';
  strcpy(phvelname,root);
  strcat(phvelname,"SAC_PHP");
  */
      if(n == 0 || n != 13) break;
      j = j + 1;
  }
  nin = j;
  fclose(in);
 
  #pragma omp parallel
  {

  #pragma omp for private (ii)
  for (ii=0;ii<nin;ii++) {
   fprintf(stderr,"%dth !!!\n",ii);
   int n, npoints, nfin, nfout1, nfout2, ierr, nprpv;
   double t0, dt, delta, vmin, vmax, tmin, tmax;
   double snr, tresh, ffact, perc, taperl,fmatch,piover4;
   float sei[32768],sei_in[32768];
   double arr1[100][8],arr2[100][7];
   double c_per[100],g_vel[100],amp_p[100],snr_p[100],snr_pc[100];
   double tamp, ampo[100][32768], pred[2][300];
   int nrow, ncol, npred;
   double prpvper[300],prpvvel[300]; /* phase vel prediction files  */

   char  *p,name[160],buf[200],phvelname[300],root[160],amp_name[160];
   FILE  *inv, *fas;
   int   i,pflag;
   int   sac = 1; /* =1 - SAC, =0 - ftat files    */


   piover4 = piover4s[ii];
   vmin = vmins[ii];
   vmax = vmaxs[ii];
   tmin = tmins[ii];
   tmax = tmaxs[ii];
   tresh = treshs[ii];
   ffact = ffacts[ii];
   taperl = taperls[ii];
   snr = snrs[ii];
   fmatch = fmatchs[ii];
   pflag = pflags[ii];
   strcpy(name,names[ii]);
   strcpy(phvelname,phvelnames[ii]);

//      printf("pi/4= %4.1lf, vmin= %lf, vmax= %lf, tmin= %lf, tmax= %lf\n",
//              piover4,vmin,vmax,tmin,tmax);
//      printf("Tresh= %lf, Filter factor= %lf, taperl= %lf, SNR= %lf, Match = %lf\nData file name=%s\n",
//             tresh,ffact,taperl,snr,fmatch,name);
/* if presents, read phase velocity prediction file
   ---  */
  nprpv = 0;
  if((inv = fopen(phvelname,"r")) == NULL) {
      printf("Can not find file %s.\n",phvelname);
  } else {
  while(fgets(buf,200,inv) != NULL) {
         if((n = sscanf(buf,"%lf %lf",&prpvper[nprpv],&prpvvel[nprpv])) < 2) break;
         nprpv++;
     }
         fclose(inv);
//         printf("Phase velocity prediction file name= %s\n",phvelname);
  }

/*
 *   read SAC or ascii data 
 */
      readdata_(&sac,name,&n,&dt,&delta,&t0,sei,&ierr);
  nfin    = 64;
  npoints = 5;          /* only 3 points in jump              */
  perc    = 50.0;       /* 50 % for output segment            */
//  printf("#filters= %d, Perc= %6.2f %s, npoints= %d, t0= %6.2f\n",
//          nfin,perc,"%",npoints,t0);

//  printf("FTAN - the first ineration\n");
/*  nprpv = 0;    */
  ffact =1.0;

// put data in sei_in //
  //fprintf(stderr,"%d \n",n);
  //abort();
  int flag = n;
  for (i=0;i<n;i++) { sei_in[i] = sei[i]; if (sei[i]==0.) flag = flag - 1; }
  if (flag < n/10.) { fprintf(stderr,"lot's of zeros! %s\n",name);continue;  }
///////////////////
  fprintf(stderr,"now do aftan1 %d %dth thread\n",ii,omp_get_thread_num());
  aftanpg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh,
        &ffact,&perc,&npoints,&taperl,&nfin,&snr,&nprpv,prpvper,prpvvel,
        &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);
  printres_(&dt,&delta,&nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,
            ampo,&ierr,name,"_1",pflag);
            //ampo,&ierr,name,"_1",pflag);
  if(nfout2 == 0) continue;   /* break aftan sequence     */
//  printf("Tamp = %9.3lf, nrow = %d, ncol = %d\n",tamp,nrow,ncol);

/* Make prediction based on the first iteration               */

  npred = nfout2;
  tmin = arr2[0][1];
  tmax = arr2[nfout2-1][1];
  for(i = 0; i < nfout2; i++) {
      pred[0][i] = arr2[i][1];   /* apparent period  */
      pred[1][i] = arr2[i][2];   /* group velocities */
  }

/* FTAN with phase with phase match filter. Second Iteration. */
  printf("FTAN - the second iteration (phase match filter)\n");
  ffact = 1.0;
//  printf("Filter factor=%6.2lf\n",ffact);
  aftanipg_(&piover4,&n,sei,&t0,&dt,&delta,&vmin,&vmax,&tmin,&tmax,&tresh,
        &ffact,&perc,&npoints,&taperl,&nfin,&snr,&fmatch,&npred,pred,
        &nprpv,prpvper,prpvvel,
        &nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,ampo,&ierr);
  //for (i=0;i<nfout2;i++)fprintf(stderr,"%g %g\n",arr2[i][2],arr2[i][3]);
  printres_(&dt,&delta,&nfout1,arr1,&nfout2,arr2,&tamp,&nrow,&ncol,
            ampo,&ierr,name,"_2",pflag);
            //ampo,&ierr,name,"_2",pflag);
//  printf("Tamp = %9.3lf, nrow = %d, ncol = %d\n",tamp,nrow,ncol);

  //fprintf(stderr,"nnn %d\n",n);
  //}

  /////////////////////////////// snr and amplitude /////////////////////////////////////
  for(i = 0; i < nfout2; i++) {
      c_per[i]=arr2[i][0];
      g_vel[i]=arr2[i][2];
     //fprintf(stderr,"%g %g\n",c_per[i],arr2[i][3]);
     }
  
  //fprintf(stderr,"%d\n",n);
  //abort();
  get_snr(sei_in,n,dt,delta,t0,c_per,g_vel,nfout2,amp_p,snr_p,snr_pc); // positive lag
  //fprintf(stderr,"after get_snr: %g %s\n",amp_p[0],name);
  sprintf(amp_name,"%s_amp_snr",name);
  if((fas=fopen(amp_name,"w"))==NULL) {
     printf("Cannot open file %s to write!\n", amp_name);
     exit (1);
    }
  for(i = 0; i < nfout2; i++)
     fprintf(fas,"%8.4f   %.5g  %8.4f %8.4f\n",arr2[i][1],amp_p[i],snr_p[i],snr_pc[i]);
  fclose(fas);
  ///////////////////////////////////////////////////////////////////////////////////
  }
 }
  //fclose(in);
  return 0;
}
