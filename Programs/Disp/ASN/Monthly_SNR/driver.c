#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "mysac64.h"
//#include "/home/tianye/MyLib/Dis_Azi.h"
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_num_threads() { return 1;}
inline omp_int_t omp_set_dynamic(int) {return 0;}
inline omp_int_t omp_set_num_threads(int) {return 0;}
#endif

#define NSTA 2000
#define BLKs 5000
#define NFRQ 300

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD);

int get_snr(float *sei, int nsample, double dt, double dist, double b, double *c_per, double *g_vel, int nper, double *amp_max, double *snr2, double *noise);

int ComputeSNR(char *dispname, char *outname, float * sei_p, SAC_HD *shd) {
   /* check/open dispersion file for positive lag */
   char buff[300];
   int i, nper = 0;
   double ftmp, cper[NFRQ], gvel[NFRQ], amp[NFRQ], noise[NFRQ], snr[NFRQ], phvtmp;
   FILE *fdisp;
   if( (fdisp=fopen(dispname, "r")) == NULL ) {
      std::cerr<<"Error(ComputeSNR): Cannot access dispersion file "<<dispname<<std::endl;
      return 0;
   }
   else { /* read in dispersion information */
      for(i=0; fgets(buff, 300, fdisp); i++ ) {
         if( sscanf(buff, "%lf %lf %lf %lf %lf", &ftmp, &ftmp, &cper[i], &gvel[i], &phvtmp) != 5 ) continue;
      }
      fclose(fdisp);
      nper = i;
   }

   /* compute snr on the positive lag */
   int n = shd->npts; 
   double dt = shd->delta, dist = shd->dist, t0 = shd->b;
   get_snr(sei_p, n, dt, dist, t0, cper, gvel, nper, amp, snr, noise); // positive lag

   /* output amp and snr for the positive lag */
   FILE * fout;
   if((fout=fopen(outname,"w"))==NULL) {
      printf("Error(main): Cannot write to file %s\n", outname);
      exit (0);
   }
   for(i = 0; i < nper; i++) fprintf(fout,"%8.4f   %.5g  %8.4f  %.5g  %.0f\n", cper[i], amp[i], snr[i], noise[i], shd->user0);
   fclose(fout);

   return 1;
}


int main(int argc, char *argv[]) 
{
   /* check input parameters */
   if( argc != 4 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [sac file (monthly CC)] [Disp file pos (from FTAN)] [Disp file neg (from FTAN)]"<<std::endl;
      exit(-1);
   }

   /* read sac signal and header */
   SAC_HD shd;
   float *sei_p = NULL;
   if( read_sac(argv[1], &sei_p, &shd) == NULL ) {
      fprintf(stderr, "ERROR(read_sac): %s\n", argv[1]);
      exit(-1);
   }
   if( shd.b>0. || shd.e<0. ) {
      std::cerr<<"Error(main): Expecting symetric sac record!"<<std::endl;
      exit(-1);
   }
   /* and split into sei_p and sei_n */
   int i, len = shd.npts/2;
   int n = len+1;
   float sei_n[n];
   for(i=0;i<n;i++) sei_n[i]=sei_p[len-i];
   for(i=0;i<n;i++) sei_p[i]=sei_p[len+i];
   //for(k=0;k<=len;k++) sei[k]=(sei_p[k]+sei_n[k])/2.;
   shd.npts = n; shd.b = 0.;

   /* compute SNR for positive lag */
   char outname[300];
   sprintf(outname, "%s_pos_amp_snr", argv[1]);
   ComputeSNR( argv[2], outname, sei_p, &shd );

   /* compute SNR for negative lag */
   sprintf(outname, "%s_neg_amp_snr", argv[1]);
   ComputeSNR( argv[2], outname, sei_n, &shd );

   free(sei_p);
   return 0;
}
