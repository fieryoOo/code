/*
This code extracts information (distance, azimuth, velocity, amplitude, and snr of both positive and negative lags) from FTAN results.
   a measurement is accepted if it matches the input predicted dispersions 
Input 1: station.lst (STANM LON LAT)
Input 2: input-file list ( staname1 staname2 fDisppos fAmppos fDispneg fAmpneg daynum fpredDispgrv fpredDispphv )
*/

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include "/home/tianye/MyLib/Dis_Azi.h"
using namespace std;

#define NSTA 1000
#define NSPR 100000
#define NFRQ 100

struct STATION {
   char name[6];
   float lon, lat;
};

struct STAPAIR {
   char sta1[6], sta2[6];
   char disp_pf[NFRQ], disp_nf[NFRQ];
   char snr_pf[NFRQ], snr_nf[NFRQ];
   char pdisp_g[NFRQ], pdisp_p[NFRQ];
   float daynum;
};

struct DATA2C {
   float x[NFRQ];
   float y[NFRQ];
   int ndat;
};

/* global data array */
struct DATA2C grvpos, phvpos, amppos, snrpos;
struct DATA2C grvneg, phvneg, ampneg, snrneg;
struct DATA2C grvpre, phvpre;

void InserSort(struct DATA2C *dat) {
   int i, j;
   float ftmpx, ftmpy;
   for(i=1;i<dat->ndat;i++) {
      ftmpx = dat->x[i]; ftmpy = dat->y[i];
      for(j=i; j>0 && ftmpx<dat->x[j-1]; j--) { dat->x[j] = dat->x[j-1]; dat->y[j] = dat->y[j-1]; }
      dat->x[j] = ftmpx; dat->y[j] = ftmpy;
   }
}

void InserSort1(float *arr, float *dat1, int n) {
   int i,j;
   float ftmp0, ftmp1;
   for(i=1;i<n;i++) {
      ftmp0=arr[i]; ftmp1=dat1[i];
      for(j=i;j>0 && ftmp0<arr[j-1]; j--) { arr[j]=arr[j-1]; dat1[j]=dat1[j-1]; }
      arr[j]=ftmp0; dat1[j]=ftmp1;
   }
}

void InserSort2(float *arr, float *dat1, float *dat2, int n) {
   int i,j;
   float ftmp0, ftmp1, ftmp2;
   for(i=1;i<n;i++) {
      ftmp0=arr[i]; ftmp1=dat1[i]; ftmp2=dat2[i];
      for(j=i;j>0 && ftmp0<arr[j-1]; j--) { arr[j]=arr[j-1]; dat1[j]=dat1[j-1]; dat2[j]=dat2[j-1]; }
      arr[j]=ftmp0; dat1[j]=ftmp1; dat2[j]=ftmp2;
   }
}

int readdispO( char * fdisp, char *fsnr, struct DATA2C *grv, struct DATA2C *phv, struct DATA2C *amp, struct DATA2C *snr ) {
   FILE *fin = NULL;
   char buff[300];
   int i, n;
   float ftmp;
   // read in and sort observed velocities
   if((fin=fopen(fdisp,"r"))==NULL) return 0;
   for(i=0;fgets(buff, 300, fin)!=NULL;) {
      if(sscanf(buff, "%f %f %f %f %f", &ftmp, &ftmp, &(grv->x[i]), &(grv->y[i]), &(phv->y[i]))!=5) continue;
      phv->x[i] = grv->x[i++];
   }
   fclose(fin);
   grv->ndat = i; phv->ndat = i;
   InserSort(grv); InserSort(phv);
}
//      readdispP( spr[ipth].pdisp_g, spr[ipth].pdisp_p, &grvpre, &phvpre );

int get_amp_snr(char *fname, float pero, float *ampo, float *snro) {
   FILE *ff;
   char buff[300];
   int i, n;
   float per[NFRQ], amp[NFRQ], snr[NFRQ];
   *amp = -1, *snro = -1;
   if((ff=fopen(fname,"r"))==NULL) return 0;
   for(i=0;fgets(buff, 300, ff)!=NULL;i++)
      if(sscanf(buff, "%f %f %f", &per[i], &amp[i], &snr[i])!=3) continue;
   fclose(ff);
   n = i;
   InserSort2(per,amp,snr,n);
   for(i=0;i<n && per[i]<pero; i++) {}
   if(i==0 || i==n) return 0;
   *ampo = amp[i-1]+(pero-per[i-1])/(per[i]-per[i-1])*(amp[i]-amp[i-1]);
   *snro = snr[i-1]+(pero-per[i-1])/(per[i]-per[i-1])*(snr[i]-snr[i-1]);
   return 1;
}

int get_vel(char *fname, float pero, float *grvo, float *phvo) {
// for calculate group and phase vel from FTAN-produced DISP file
// suppose per to be at column 3 and group/phase vels at column 4/5
   FILE *ff;
   char buff[300];
   int i, n, itmp;
   float per[NFRQ], grv[NFRQ], phv[NFRQ], ftmp;
   *grvo = -1; *phvo = -1;
   if((ff=fopen(fname,"r"))==NULL) return 0;
   for(i=0;fgets(buff, 300, ff)!=NULL;i++)
      if(sscanf(buff, "%d %f %f %f %f", &itmp, &ftmp, &per[i], &grv[i], &phv[i])!=5) continue;
   fclose(ff);
   n = i;
   InserSort2(per,grv,phv,n);
   for(i=0;i<n && per[i]<pero; i++) {}
   if(i==0 || i==n) return 0;
   *grvo = grv[i-1]+(pero-per[i-1])/(per[i]-per[i-1])*(grv[i]-grv[i-1]);
   *phvo = phv[i-1]+(pero-per[i-1])/(per[i]-per[i-1])*(phv[i]-phv[i-1]);
   return 1;
}

int get_vel_pred(char *fname, float pero, float *velo) {
// for calculate vel from predicted DISP file
// suppose per to be at column 1 and vels at column 2
   FILE *ff;
   char buff[300];
   int i, n;
   float per[NFRQ], vel[NFRQ];
   *velo = -1;
   if((ff=fopen(fname,"r"))==NULL) return 0;
   for(i=0;fgets(buff, 300, ff)!=NULL;i++)
      if(sscanf(buff, "%f %f", &per[i], &vel[i])!=2) continue;
   fclose(ff);
   n = i;
   InserSort1(per,vel,n);
   for(i=0;i<n && per[i]<pero; i++) {}
   if(i==0 || i==n) return 0;
   *velo = vel[i-1]+(pero-per[i-1])/(per[i]-per[i-1])*(vel[i]-vel[i-1]);
   return 1;
}

int disp_match(char *dispname, char *grvname, char *phvname, float pero, float grvmisf, float phvmisf) { 
/* --------- result -----------
-1: can't tell
 0: do not match
 1: groups do match, phase doesn't
 2: both match
------------------------------*/
   int i, ncheck = 0, result = 2;
   float per[5];
   float grvo[5], phvo[5], grvp[5], phvp[5];
   per[0] = pero/1.1; per[1] = pero/1.05; per[2] = pero;
   per[3] = pero*1.05; per[4] = pero*1.1;
   for(i=0;i<5;i++) {
      if( !get_vel(dispname, per[i], &grvo[i], &phvo[i]) ) continue;
      if( !get_vel_pred(grvname, per[i], &grvp[i]) ) continue;
      if( !get_vel_pred(phvname, per[i], &phvp[i]) ) continue;
      //cout<<i<<": "<<grvo[i]<<" = "<<grvp[i]<<" ?  "<<phvo[i]<<" = "<<phvp[i]<<" ?"<<endl;
      if(fabs(grvo[i]-grvp[i])/grvp[i]>grvmisf) return 0;
      if(fabs(phvo[i]-phvp[i])/phvp[i]>phvmisf) result = 1;
      ncheck++;
   }
   if( ncheck < 3 ) return -1;
   return result;
}

int DistPointCurve() {
   
}
/*
int DispQuality(char *dispname, char *grvname, char *phvname, float grvmisf, float phvmisf) {
// read in group and phase disp data
   FILE *fin;
   char buff[300];
   int i, n;
   float per[NFRQ], grvel[NFRQ], phvel[NFRQ], ftmp;
   if((ff=fopen(dispname,"r"))==NULL) return 0;
   for(i=0;fgets(buff, 300, ff)!=NULL;i++)
      if(sscanf(buff, "%f %f %f %f %f", &ftmp, &ftmp, &per[i], &grvel[i], &phvel[i])!=5) continue;
   fclose(ff);
   n = i;
   InserSort2(per,grvel,phvel,n);
// for each data point, measure minimum distance to predicted curves
   DistPoint2Curve();
// throw out points with misfit larger than max-allowed misfit

// analysis remaining points and compute a score

   exit(-2);
}
*/


main(int na, char *arg[])
{
   if(na!=6) {
      cout<<"usage: "<<arg[0]<<" [station.lst] [stanm1-stanm2-fDisppos-fAmppos-fDispneg-fAmpneg-daynum-fpredDispgrv-fpredDispphv list] [period] [max grvmisf perc] [max phvmisf perc]"<<endl;
      return 0;
   }
 
   FILE *ff;
   char buff[300];
   int i, nsta, npth;
   struct STATION sta[NSTA];
   struct STAPAIR spr[NSPR];
/* read in sta.name, sta.lon, sta.lat from station list */
   if((ff=fopen(arg[1],"r"))==NULL) {
      cout<<"Cannot open file "<<arg[1]<<endl;
      return -1;
   }
   for(i=0;fgets(buff, 300, ff)!=NULL;i++) {
      sscanf(buff,"%s %f %f", sta[i].name, &sta[i].lon, &sta[i].lat);
      if(sta[i].lon<0) sta[i].lon += 360.;
   }
   nsta=i;
   fclose(ff);

/* read in sta1, sta2, dispf, snrf from file list */
   if((ff=fopen(arg[2],"r"))==NULL) {
      cout<<"Cannot open file "<<arg[2]<<endl;
      return -1;
   }
   for(i=0;fgets(buff, 300, ff)!=NULL;i++)
      sscanf(buff,"%s %s %s %s %s %s %f %s %s", spr[i].sta1, spr[i].sta2, spr[i].disp_pf, spr[i].snr_pf, spr[i].disp_nf, spr[i].snr_nf, &(spr[i].daynum), spr[i].pdisp_g, spr[i].pdisp_p);
   npth = i;
   fclose(ff);
/* main loop. Process data from each path. Compute/check for dist and snr, analyse dispersion quality and construct maps at each period */
   int ipth, isp, isn;
   int flagp, flagn;
   double dist, azi1, azi2;
   float grvp, grvn, phvp, phvn, ampp, ampn, snrp, snrn;
   float pero = atof(arg[3]), grvmisf = atof(arg[4]), phvmisf = atof(arg[5]);

   sprintf(buff, "Disp_info_wp_%.1fsec", pero);
   ff = fopen(buff, "w");
   fprintf(ff, "sta1 lat1 lon1  sta2 lat2 lon2  dist azi1 azi2 dnum : snrp ampp/dnum grvp phvp  snrn ampn/dnum grvn phvn\n");
   for(ipth=0;ipth<npth;ipth++) {
   // Search for sta1 and sta2 in station list
      for(isp=0;isp<nsta;isp++) if(strcmp(spr[ipth].sta1, sta[isp].name)==0) break;
      if(isp==nsta) continue;
      for(isn=0;isn<nsta;isn++) if(strcmp(spr[ipth].sta2, sta[isn].name)==0) break;
      if(isn==nsta) continue;
      cout<<"Extracting information for path "<<spr[ipth].sta1<<" - "<<spr[ipth].sta2<<endl;
      calc_dist(sta[isp].lat, sta[isp].lon, sta[isn].lat, sta[isn].lon, &dist);
      calc_azimuth(sta[isp].lat, sta[isp].lon, sta[isn].lat, sta[isn].lon, &azi1);
      calc_azimuth(sta[isn].lat, sta[isn].lon, sta[isp].lat, sta[isp].lon, &azi2);
   // read in data ( observed grv, phv, amp, snr for pos and neg component; predicted grv, phv )
      //readdispO( spr[ipth].disp_pf, spr[ipth].snr_pf, &grvpos, &phvpos, &amppos, &snrpos );
      //readdispO( spr[ipth].disp_nf, spr[ipth].snr_nf, &grvneg, &phvneg, &amppos, &snrneg );
      //readdispP( spr[ipth].pdisp_g, spr[ipth].pdisp_p, &grvpre, &phvpre );
   // Check for qualitiy of each measured dispersion by comparing it with the predicted curves.
    //Positive lag:
      //DispQuality(spr[ipth].disp_pf, spr[ipth].pdisp_g, spr[ipth].pdisp_p, grvmisf, phvmisf);
      flagp = disp_match(spr[ipth].disp_pf, spr[ipth].pdisp_g, spr[ipth].pdisp_p, pero, grvmisf, phvmisf);
      if( flagp<=0 ) {
	 ampp = -1; snrp = -1;
	 grvp = -1; phvp = -1;
      }
      else {
	 get_amp_snr(spr[ipth].snr_pf, pero, &ampp, &snrp); //get amplitude and SNR
	 get_vel(spr[ipth].disp_pf, pero, &grvp, &phvp); //get group and phase vel
	 if( flagp==1 ) phvp = -1;
      }
    //Negative lag:
      flagn = disp_match(spr[ipth].disp_nf, spr[ipth].pdisp_g, spr[ipth].pdisp_p, pero, grvmisf, phvmisf);
      if( flagn<=0 ) {
         ampn = -1; snrn = -1;
         grvn = -1; phvn = -1;
      }
      else {
         get_amp_snr(spr[ipth].snr_nf, pero, &ampn, &snrn); //get amplitude and SNR
         get_vel(spr[ipth].disp_nf, pero, &grvn, &phvn); //get group and phase vel
         if( flagn==1 ) phvn = -1;
      }
      fprintf(ff, "%s %f %f  %s %f %f  %lf %lf %lf %f : %f %g %f %f  %f %g %f %f\n", sta[isp].name, sta[isp].lat, sta[isp].lon, sta[isn].name, sta[isn].lat, sta[isn].lon, dist, azi1, azi2, spr[ipth].daynum, snrp, ampp/spr[ipth].daynum, grvp, phvp, snrn, ampn/spr[ipth].daynum, grvn, phvn);
   }
   fclose(ff);
  
   return 1;
}
