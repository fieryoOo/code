/*
This code extracts information (distance, azimuth, velocity, amplitude, and snr of both positive and negative lags) from FTAN results.
   In general, the lag(s) with SNR lower than snrmin are considered bad measurement(s) and are rejected.
   In the case where both lags has higher than standard SNR, the lower-SNR lag is rejected if a phase mis-matching is detected
Input 1: station.lst (STANM LON LAT)
Input 2: input-file list (STA1 STA2 Disp_File_Positive SNR_File_Positive Disp_File_Negative SNR_File_Negative)
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
#define BLKSIZE 1000

struct STATION {
   char name[6];
   float lon, lat;
};

struct STAPAIR {
   char sta1[6], sta2[6];
   char disp_pf[100], disp_nf[100];
   char snr_pf[100], snr_nf[100];
   float daynum;
};

void InserSort2(float *arr, float *dat1, float *dat2, int n) {
   int i,j;
   float ftmp0, ftmp1, ftmp2;
   for(i=1;i<n;i++) {
      ftmp0=arr[i]; ftmp1=dat1[i]; ftmp2=dat2[i];
      for(j=i;j>0 && ftmp0<arr[j-1]; j--) { arr[j]=arr[j-1]; dat1[j]=dat1[j-1]; dat2[j]=dat2[j-1]; }
      arr[j]=ftmp0; dat1[j]=ftmp1; dat2[j]=ftmp2;
   }
}

int get_amp_snr(char *fname, float pero, float *ampo, float *snro) {
   FILE *ff;
   char buff[300];
   int i, n;
   float per[100], amp[100], snr[100];
   *snro = 0;
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
   FILE *ff;
   char buff[300];
   int i, n, itmp;
   float per[100], grv[100], phv[100], ftmp;
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

int ph_match(char *nfname, char *pfname, float pero, float vmisf) {
   int i;
   float per[3];
   float grvn[3], phvn[3], grvp[3], phvp[3];
   per[0] = pero/1.05, per[1] = pero, per[2] = pero/0.95;
   for(i=0;i<3;i++) {
      get_vel(nfname, per[i], &grvn[i], &phvn[i]);
      get_vel(pfname, per[i], &grvp[i], &phvp[i]);
      if(fabs(phvn[i]-phvp[i])/phvn[i]>vmisf) return 0;
   }
   return 1;
}


main(int na, char *arg[])
{
   if(na!=6) {
      cout<<"usage: "<<arg[0]<<" [station.lst] [sta1-sta2-Dispposf-SNRp-Dispnegf-SNRn-daynum list] [period] [min snr] [max vmisf perc]"<<endl;
      return 0;
   }
 
   FILE *ff;
   char buff[500];
   int i, nsta, npth;
   struct STATION sta[NSTA];
   struct STAPAIR *spr = NULL;
//read in sta.name, sta.lon, sta.lat from station list
   if((ff=fopen(arg[1],"r"))==NULL) {
      cout<<"Cannot open file "<<arg[1]<<endl;
      return -1;
   }
   for(i=0;fgets(buff, 500, ff)!=NULL;i++) {
      sscanf(buff,"%s %f %f", sta[i].name, &sta[i].lon, &sta[i].lat);
      if(sta[i].lon<0) sta[i].lon += 360.;
   }
   nsta=i;
   fclose(ff);

//read in sta1, sta2, dispf, snrf from file list
   if((ff=fopen(arg[2],"r"))==NULL) {
      cout<<"Cannot open file "<<arg[2]<<endl;
      return -1;
   }
   int nblk = 0;
   for(i=0;fgets(buff, 500, ff)!=NULL;i++) {
      if( nblk*BLKSIZE <= i ) spr = (struct STAPAIR *) realloc ( spr, (++nblk)*BLKSIZE * sizeof(struct STAPAIR));
      sscanf(buff,"%s %s %s %s %s %s %f", spr[i].sta1, spr[i].sta2, spr[i].disp_pf, spr[i].snr_pf, spr[i].disp_nf, spr[i].snr_nf, &(spr[i].daynum));
   }
   npth = i;
   fclose(ff);
//read in dispersion data from each path, check for dist and snr
   int ipth, isp, isn;
   int flagp, flagn;
   double dist, azi1, azi2;
   float grvp, grvn, phvp, phvn, ampp, ampn, snrp, snrn;
   float pero = atof(arg[3]), snrmin = atof(arg[4]), vmisf = atof(arg[5]);
   sprintf(buff, "Disp_info_%.1fsec", pero);
   ff = fopen(buff, "w");
   fprintf(ff, "sta1(1) lat1(2) lon1(3)  sta2(4) lat2(5) lon2(6)  dist(7) azi1(8) azi2(9) dnum(10) : snrp(12) ampp/dnum(13) grvp(14) phvp(15)  snrn(16) ampn/dnum(17) grvn(18) phvn(19)\n");
   for(ipth=0;ipth<npth;ipth++) {
      //search for sta1 and sta2 in station list
      for(isp=0;isp<nsta;isp++) if(strcmp(spr[ipth].sta1, sta[isp].name)==0) break;
      if(isp==nsta) continue;
      for(isn=0;isn<nsta;isn++) if(strcmp(spr[ipth].sta2, sta[isn].name)==0) break;
      if(isn==nsta) continue;
      cout<<"Extracting information for path "<<spr[ipth].sta1<<" - "<<spr[ipth].sta2<<endl;
      calc_dist(sta[isp].lat, sta[isp].lon, sta[isn].lat, sta[isn].lon, &dist);
      calc_azimuth(sta[isp].lat, sta[isp].lon, sta[isn].lat, sta[isn].lon, &azi1);
      calc_azimuth(sta[isn].lat, sta[isn].lon, sta[isp].lat, sta[isp].lon, &azi2);
      //get amplitude and SNR
      ampp = 0; snrp = 0;
      get_amp_snr(spr[ipth].snr_pf, pero, &ampp, &snrp);
      ampn = 0; snrn = 0;
      get_amp_snr(spr[ipth].snr_nf, pero, &ampn, &snrn);
      //compare pos and neg lags to decide which lag has meaningful signal
      flagp = 0; flagn = 0;
      if( snrp>snrmin && snrn>snrmin ) {
         if( ph_match(spr[ipth].disp_nf, spr[ipth].disp_pf, pero, vmisf) ) { flagp = 1; flagn = 1; }
	 else if (snrp > snrn) flagp = 1;
	 else flagn = 1;
      }
      else if (snrp > snrmin) flagp = 1;
      else if (snrn > snrmin) flagn = 1;
      //get grv and phv from the meaningful lag
      grvp = 0; phvp = 0; grvn = 0; phvn = 0;
      if( flagp ) get_vel(spr[ipth].disp_pf, pero, &grvp, &phvp); 
      else snrp = 0;
      if( flagn ) get_vel(spr[ipth].disp_nf, pero, &grvn, &phvn); 
      else snrn = 0;
      fprintf(ff, "%s %f %f  %s %f %f  %lf %lf %lf %f : %f %g %f %f  %f %g %f %f\n", sta[isp].name, sta[isp].lat, sta[isp].lon, sta[isn].name, sta[isn].lat, sta[isn].lon, dist, azi1, azi2, spr[ipth].daynum, snrp, ampp/spr[ipth].daynum, grvp, phvp, snrn, ampn/spr[ipth].daynum, grvn, phvn);
   }
   fclose(ff);
 
   free(spr); 
   return 1;
}
