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

float perl, perh, ghlgrv, ghlphv, alphagrv, alphaphv;

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

void InserSort1(float *arr, float *dat1, int n) {
   int i,j;
   float ftmp0, ftmp1;
   for(i=1;i<n;i++) {
      ftmp0=arr[i]; ftmp1=dat1[i];
      for(j=i;j>0 && ftmp0<arr[j-1]; j--) { arr[j]=arr[j-1]; dat1[j]=dat1[j-1]; }
      arr[j]=ftmp0; dat1[i]=ftmp1;
   }
}

void InserSort2(float *arr, float *dat1, float *dat2, int n) {
   int i,j;
   float ftmp0, ftmp1, ftmp2;
   for(i=1;i<n;i++) {
      ftmp0=arr[i]; ftmp1=dat1[i]; ftmp2=dat2[i];
      for(j=i;j>0 && ftmp0<arr[j-1]; j--) { arr[j]=arr[j-1]; dat1[j]=dat1[j-1]; dat2[j]=dat2[j-1]; }
      arr[j]=ftmp0; dat1[i]=ftmp1; dat2[i]=ftmp2;
   }
}

// compute distance% in the freq-vel domain
float FVDist(float per0, float vel0, float *per, float *vel, int n) {
   float dx, dy, dis, dismin = 9999999.;
   int i;
   for(i=0;i<n;i++) {
      dy = 1.-vel0/vel[i];
      dx = 1.-per0/per[i];
      dis = dx*dx + dy*dy;
      if(dismin>dis) dismin = dis;
   }
   return sqrt(dismin);
}

int AmpSNRSummation( char *fdisp, char *fsnr, char *fgrvpred, char *fphvpred, float *snr, float *amp ) {
   FILE *fin1, *fin2;
   char buff1[300], buff2[300];
   int i, nfgrv, nfphv;
   float pergrv[NFRQ], grvpred[NFRQ], perphv[NFRQ], phvpred[NFRQ];
   // read predicted group disp
   if((fin1=fopen(fgrvpred, "r")) == NULL) return 0;
   for(i=0;fgets(buff1, 300, fin1)!=NULL;) {
      if(sscanf(buff1, "%f %f", &pergrv[i], &grvpred[i])!=2) continue;
      i++;
   }
   fclose(fin1);
   nfgrv = i;
   // define period range
   float perlow = pergrv[0], perhigh = pergrv[0];
   for(i=1; i<nfgrv; i++) {
      if( pergrv[i] < perlow ) perlow = pergrv[i];
      else if (pergrv[i] > perhigh ) perhigh = pergrv[i];
   }
   if( perlow<perl ) perlow = perl;
   if( perhigh>perh ) perhigh = perh;
   // read predicted phase disp and define period range
   if((fin1=fopen(fphvpred, "r")) == NULL) return 0;
   for(i=0;fgets(buff1, 300, fin1)!=NULL;) {
      if(sscanf(buff1, "%f %f", &perphv[i], &phvpred[i])!=2) continue;
      i++;
   }
   fclose(fin1);
   nfphv = i;
   // read observed parameters and weighted-sum amp and snr
   *snr = 0.; *amp = 0.;
   int ndat = 0, neff = 0;
   float ftmp, weit = 0., weight, disgrv, disphv;
   float percur, grvcur, phvcur, snrcur, ampcur;
   float perlowO = 999999., perhighO = 0.;
   if((fin1=fopen(fdisp, "r")) == NULL || ((fin2=fopen(fsnr, "r"))==NULL)) return 0;
   while(fgets(buff1, 300, fin1) && fgets(buff2, 300, fin2)) {
      if( sscanf(buff1, "%f %f %f %f %f", &ftmp, &ftmp, &percur, &grvcur, &phvcur) != 5 ) continue;
      if( sscanf(buff2, "%f %f %f %f %f", &ftmp, &ampcur, &snrcur) != 3 ) continue;
      if( ftmp != percur ) {
	 cerr<<fdisp<<" and "<<fsnr<<" does not match!"<<endl;
	 continue;
      }
      if( percur<perlow || percur>perhigh ) continue;
      if( percur < perlowO ) perlowO = percur;
      if( percur > perhighO ) perhighO = percur;
      ndat++;
      disgrv = FVDist(percur, grvcur, pergrv, grvpred, nfgrv);
      if( disgrv > ghlgrv*3 ) continue;
      disphv = FVDist(percur, phvcur, perphv, phvpred, nfphv);
      neff++;
      weight = exp(alphagrv*disgrv*disgrv) * 0.7 + exp(alphaphv*disphv*disphv) * 0.3;
      *snr += weight*snrcur;
      *amp += weight*ampcur;
      weit += weight;
//cerr<<fdisp<<" "<<disgrv<<" "<<disphv<<" "<<*snr<<" "<<*amp<<" "<<weit<<endl;
   }
   fclose(fin1); fclose(fin2);
   // compute effective signal snr and amp in the period range
   if( neff == 0 ) return 1;
   float effperc = (log(perhighO)-log(perlowO))/(log(perhigh)-log(perlow))*neff/ndat/weit;
   *snr *= effperc; *amp *= effperc;
   
   return 1;
}

main(int na, char *arg[])
{
   if(na!=7) {
      cout<<"usage: "<<arg[0]<<" [station.lst] [stanm1-stanm2-fDisppos-fAmppos-fDispneg-fAmpneg-daynum-fpredDispgrv-fpredDispphv list] [period_lowerend] [period_higherend] [ghlgrv] [ghlphv]"<<endl;
      return 0;
   }
 
/* ghlgrv & ghlphv
Gaussian halflength as percentage for group and phase dispersion shift
defines how the amplitude and SNRs are weighted
*/
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
/* main loop. Process data from each path. Compute dist and snr, and as an indicator of FTAN quality,
   compute summation of SNR&amp weighted by distance between observed and predicted dispersion curves */
   int ipth, isp, isn;
   int flagp, flagn;
   float snrsig_pos, snrsig_neg, ampsig_pos, ampsig_neg;
   double dist, azi1, azi2;
   perl = atof(arg[3]); perh = atof(arg[4]); 
   ghlgrv = atof(arg[5]); ghlphv = atof(arg[6]);
   alphagrv = -0.5/(ghlgrv*ghlgrv); alphaphv = -0.5/(ghlphv*ghlphv);

   sprintf(buff, "Disp_quality_%f_%f", perl, perh);
   ff = fopen(buff, "w");
   fprintf(ff, "sta1 lat1 lon1  sta2 lat2 lon2  dist azi1 azi2 dnum : snrsig_pos ampsig_pos/dnum snrsig_neg ampsig_neg/dnum\n");
   for(ipth=0;ipth<npth;ipth++) {
   // Search for sta1 and sta2 in station list
      for(isp=0;isp<nsta;isp++) if(strcmp(spr[ipth].sta1, sta[isp].name)==0) break;
      if(isp==nsta) continue;
      for(isn=0;isn<nsta;isn++) if(strcmp(spr[ipth].sta2, sta[isn].name)==0) break;
      if(isn==nsta) continue;
      cout<<"Extracting information for path "<<spr[ipth].sta1<<" - "<<spr[ipth].sta2<<endl;
   // compute distance and azimuth
      calc_dist(sta[isp].lat, sta[isp].lon, sta[isn].lat, sta[isn].lon, &dist);
      calc_azimuth(sta[isp].lat, sta[isp].lon, sta[isn].lat, sta[isn].lon, &azi1);
      calc_azimuth(sta[isn].lat, sta[isn].lon, sta[isp].lat, sta[isp].lon, &azi2);
   // compute signal SNR&amp
      if( ! AmpSNRSummation( spr[ipth].disp_pf, spr[ipth].snr_pf, spr[ipth].pdisp_g, spr[ipth].pdisp_p, &snrsig_pos, &ampsig_pos ) ) continue;
      if( ! AmpSNRSummation( spr[ipth].disp_nf, spr[ipth].snr_nf, spr[ipth].pdisp_g, spr[ipth].pdisp_p, &snrsig_neg, &ampsig_neg ) ) continue;
      fprintf(ff, "%s %f %f  %s %f %f  %lf %lf %lf %f : %f %g  %f %g\n", sta[isp].name, sta[isp].lat, sta[isp].lon, sta[isn].name, sta[isn].lat, sta[isn].lon, dist, azi1, azi2, spr[ipth].daynum, snrsig_pos, ampsig_pos/spr[ipth].daynum, snrsig_neg, ampsig_neg/spr[ipth].daynum);
   }
   fclose(ff);
  
   return 1;
}
