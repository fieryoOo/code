/*
This code investigate the quality of FTAN results by computing
  1. a quality factor defined as a weighted-average of the SNR|amplitude in a given period band, and
  2. the weighted-averaged timeshift between the positive and the negative lags of each path. 
  Weights are defined based on the distance between measured and predicted dispersion curves.
  Output: (distance, azimuth, integrated_amp, and integrated_snr of positive and negative lags and
    the timeshift between the 2 lags)
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
#include "omp.h"
using namespace std;

#define NSTA 2000
#define BLKs 5000
#define NFRQ 300

float perl, perh, ghlgrv, ghlphv, alphagrv, alphaphv, snrmin;
omp_lock_t readlock;

struct STATION {
   char name[6];
   float lon, lat;
};

struct STAPAIR {
   char sta1[6], sta2[6];
   char disp_pf[150], disp_nf[150];
   char snr_pf[150], snr_nf[150];
   char pdisp_g[150], pdisp_p[150];
   float daynum;
};

struct LAG_DATA {
   float per[NFRQ], snr[NFRQ];
   float disgrv[NFRQ], disphv[NFRQ];
   float grv[NFRQ], phv[NFRQ];
   float effperc; // percentage of the effective bandwidth 
   int ndat;
};

struct TIMING_DATA {
   float logper[NFRQ*2], snreff[NFRQ*2]; // snreff is simply the smaller of the 2 lags
   float weigrv[NFRQ*2], weiphv[NFRQ*2];
   float tdiffgrv[NFRQ*2], tdiffphv[NFRQ*2];
   float Qfactor; // overall quality factor for measurements on this particular path ('Qfactor=1' == perfect)
   int npair;
};


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

// compute distance% in the freq-vel domain
float FVDist(float per0, float vel0, float *per, float *vel, int n) {
   float dx, dy, dis, logper0 = log(per0);
   int i;
   // interpolate to get the vertical distance
   /*
   for(i=0;per[i]<per0;i++){};
   if( i < n ) {
      float veli = vel[i-1] + (vel[i]-vel[i-1]) * (logper0-log(per[i-1]))/(log(per[i])-log(per[i-1]));
      dismin = fabs(veli-vel0)/veli;
   }
   else { dismin = 1.e25; }
   */
   float veli, dismin = 1.e25;
   for(i=1;i<n;i++) {
      if( per[i-1]>per0 || per[i]<per0 ) continue;
      veli = vel[i-1] + (vel[i]-vel[i-1]) * (logper0-log(per[i-1]))/(log(per[i])-log(per[i-1]));
      dis = fabs(veli-vel0)/veli;
      if( dismin > dis ) dismin = dis;
   }
   // check for smaller (non-interpolated) distance
   for(i=0;i<n;i++) {
      dx = 1.-logper0/log(per[i]);
      if( fabs(dx) > dismin ) continue;
      dy = 1.-vel0/vel[i];
      dis = dx*dx + dy*dy;
      if(dismin>dis) dismin = dis;
   }
   return sqrt(dismin);
}

int AmpSNRSummation( char *fdisp, char *fsnr, char *fgrvpred, char *fphvpred, float *snr, float *amp, struct LAG_DATA *data ) {
   FILE *fin1, *fin2;
   char buff1[300], buff2[300];
   int i, nfgrv, nfphv;
   float pergrv[NFRQ], grvpred[NFRQ], perphv[NFRQ], phvpred[NFRQ];
  // omp_set_lock(&readlock);
   // read predicted group disp
   if((fin1=fopen(fgrvpred, "r")) == NULL) return -1;
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
   if((fin1=fopen(fphvpred, "r")) == NULL) return -1;
   for(i=0;fgets(buff1, 300, fin1)!=NULL;) {
      if(sscanf(buff1, "%f %f", &perphv[i], &phvpred[i])!=2) continue;
      i++;
   }
   fclose(fin1);
   nfphv = i;
  // omp_unset_lock(&readlock);
   // read observed parameters and weighted-sum amp and snr
   *snr = 0.; *amp = 0.;
   int ndat = 0, neff = 0;
   float ftmp, weit = 0., weight, disgrv, disphv;
   float percur, grvcur, phvcur, snrcur, ampcur;
   float perlowO = 999999., perhighO = 0.;
   if((fin1=fopen(fdisp, "r")) == NULL || ((fin2=fopen(fsnr, "r"))==NULL)) return -1;
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
      disgrv = FVDist(percur, grvcur, pergrv, grvpred, nfgrv);
      disphv = FVDist(percur, phvcur, perphv, phvpred, nfphv);
      // save data
      data->per[ndat] = percur; data->snr[ndat] = snrcur;
      data->grv[ndat] = grvcur; data->phv[ndat] = phvcur;
      data->disgrv[ndat] = disgrv; data->disphv[ndat] = disphv;
      ndat++;
      if( disgrv > ghlgrv*3 || snrcur<snrmin ) continue;
      // integrate
      weight = exp(alphagrv*disgrv*disgrv)*0.7 + exp(alphaphv*disphv*disphv) * 0.3;
      *snr += weight*snrcur;
      *amp += weight*ampcur;
      weit += weight;
      neff++;
      //cerr<<fdisp<<" "<<disgrv<<" "<<disphv<<" "<<*snr<<" "<<*amp<<" "<<weit<<endl;
   }
   fclose(fin1); fclose(fin2);
   data->ndat = ndat;
   // compute effective signal snr and amp in the period range
   if( neff == 0 ) return 0;
   data->effperc = (log(perhighO)-log(perlowO))/(log(perhigh)-log(perlow))*neff/ndat;
   ftmp = data->effperc/weit;
   *snr *= ftmp; *amp *= ftmp;
   
   return 1;
}

void ComputeDiff( struct LAG_DATA *datacurve, int ilow, int ihigh, struct LAG_DATA *datapoint, int ipoint, double dist, struct TIMING_DATA *datadiff ) {
   // avoid duplications
   int idiff = datadiff->npair;
   float logper = log(datapoint->per[ipoint]), snreff = datapoint->snr[ipoint];
   if( logper == datadiff->logper[idiff-1] ) return;
   // interpolate for tdiffgrv, tdiffphv, disgrv and disphv on datacurve at per
   float logperl = log(datacurve->per[ilow]), logperh = log(datacurve->per[ihigh]);
   float frac = (logper-logperl) / (logperh-logperl);
   float snr_curve = datacurve->snr[ilow] + (datacurve->snr[ihigh] - datacurve->snr[ilow]) * frac;
   if( snr_curve < snreff ) snreff = snr_curve;
   if( snreff < snrmin ) return;
   float grv_curve = datacurve->grv[ilow] + (datacurve->grv[ihigh] - datacurve->grv[ilow]) * frac;
   float phv_curve = datacurve->phv[ilow] + (datacurve->phv[ihigh] - datacurve->phv[ilow]) * frac;
   float disgrv_curve = datacurve->disgrv[ilow] + (datacurve->disgrv[ihigh] - datacurve->disgrv[ilow]) * frac;
   float disphv_curve = datacurve->disphv[ilow] + (datacurve->disphv[ihigh] - datacurve->disphv[ilow]) * frac;
   // compute and store tdiff & weight
   float disgrv = disgrv_curve + datapoint->disgrv[ipoint];
   float disphv = disphv_curve + datapoint->disphv[ipoint];
   float ftmp = sqrt(snreff/snrmin);
   datadiff->logper[idiff] = logper;
   datadiff->snreff[idiff] = snreff;
   datadiff->weigrv[idiff] = exp(alphagrv*disgrv*disgrv)*ftmp;
   datadiff->weiphv[idiff] = exp(alphaphv*disphv*disphv)*ftmp;
   datadiff->tdiffgrv[idiff] = dist/grv_curve - dist/datapoint->grv[ipoint];
   datadiff->tdiffphv[idiff] = dist/phv_curve - dist/datapoint->phv[ipoint];
   datadiff->npair += 1;
   //cerr<<datapoint->per[ipoint]<<" "<<grv_curve<<" "<<datapoint->grv[ipoint]<<" "<<datadiff->tdiffgrv[idiff]<<" "<<snreff<<" "<<datadiff->weigrv[idiff]<<endl;
}

void AverageTimeShift(struct TIMING_DATA *datadiff, float *TShiftG, float *sigmaG, float *TShiftP, float *sigmaP) {
   //cerr<<"Before: "<<*TShiftG<<" "<<*sigmaG<<"    "<<*TShiftP<<" "<<*sigmaP<<endl;
   int i, npair = datadiff-> npair, npef;
 // group
   float V1grv = 0., snravg = 0.;
   *TShiftG = 0.;
   for(npef=0,i=0;i<npair;i++) {
      if( (datadiff->weigrv[i]) < 1.e-10 ) continue;
      *TShiftG += (datadiff->tdiffgrv[i]) * (datadiff->weigrv[i]);
      snravg += (datadiff->snreff[i]) * (datadiff->weigrv[i]);
      V1grv += (datadiff->weigrv[i]); 
      npef++;
   }
   if( npef<3 || V1grv < 1. ) { 
      *TShiftG = -1.; *sigmaG = -1.; 
      *TShiftP = -1.; *sigmaP = -1.;
      return;
   }
   *TShiftG /= V1grv; snravg /= V1grv;
   *sigmaG = 0.;
   float V2grv=0., ftmp, weight;
   for(i=0;i<npair;i++) {
      //grv
      weight = datadiff->weigrv[i];
      if( weight < 1.e-10 ) continue;
      ftmp = (datadiff->tdiffgrv[i]) - *TShiftG;
      *sigmaG += weight*ftmp*ftmp;
      V2grv += weight*weight;
   }
   *sigmaG = sqrt(*sigmaG * V1grv) / (V1grv*V1grv-V2grv); // this is std of the mean
   // sigma decrease (, which means the computed TimeShift is more beliveable, ) when snravg is high
   datadiff->Qfactor *= 2.;
   //datadiff->Qfactor = 1.; // may exclude the effect of effperc because it has already been included in computing std of the mean
   datadiff->Qfactor *= sqrt(snravg/snrmin);
   *sigmaG /= (datadiff->Qfactor);
 // phase
   float V1phv = 0.;
   *TShiftP = 0.;
   for(npef=0,i=0;i<npair;i++) {
      if( (datadiff->weiphv[i]) < 1.e-10 ) continue;
      *TShiftP += (datadiff->tdiffphv[i]) * (datadiff->weiphv[i]);
      V1phv += (datadiff->weiphv[i]);
      npef++;
   }
   if( npef < 3 || V1phv < 1. ) { 
      *TShiftP = -1.; *sigmaP = -1.; 
      return;
   }
   *TShiftP /= V1phv; 
   *sigmaP = 0.;
   float V2phv=0.;
   for(i=0;i<npair;i++) {
      //phv
      weight = datadiff->weiphv[i];
      if( weight < 1.e-10 ) continue;
      ftmp = (datadiff->tdiffphv[i]) - *TShiftP;
      *sigmaP += weight*ftmp*ftmp;
      V2phv += weight*weight;
   }
   *sigmaP = sqrt(*sigmaP * V1phv) / (V1phv*V1phv-V2phv);
   *sigmaP /= (datadiff->Qfactor);
}

void TimeShift( struct LAG_DATA *datapos, struct LAG_DATA *dataneg, double dist, float *TShiftG, float *sigmaG, float *TShiftP, float *sigmaP ) {
   int ii, ip, in, np=datapos->ndat, nn=dataneg->ndat;
   struct TIMING_DATA datadiff;
   if(datapos->effperc < dataneg->effperc) datadiff.Qfactor = datapos->effperc;
   else datadiff.Qfactor = dataneg->effperc;
   datadiff.npair = 0;
   // compute timeshift and weight at each point
   for(ip=0,in=0; ip<np&&in<nn; ) {
      if( datapos->per[ip] < dataneg->per[in] ) {
	 for(ii=in; ii<nn && dataneg->per[ii]<=datapos->per[ip+1]; ii++) ComputeDiff(datapos, ip, ip+1, dataneg, ii, dist, &datadiff);
	 if( ii == nn ) break;
	 ip++;
      }
      else {
	 for(ii=ip; ii<np && datapos->per[ii]<dataneg->per[in+1]; ii++) ComputeDiff(dataneg, in, in+1, datapos, ii, -dist, &datadiff);
	 if( ii == np ) break;
         in++;
      }
   }
   // compute averaged timeshift and uncertainty for this sta-pair
   AverageTimeShift(&datadiff, TShiftG, sigmaG, TShiftP, sigmaP);
}

main(int na, char *arg[])
{
   if(na!=8) {
      cout<<"usage: "<<arg[0]<<" [station.lst] [stanm1-stanm2-fDisppos-fAmppos-fDispneg-fAmpneg-daynum-fpredDispgrv-fpredDispphv list] [period_lowerend] [period_higherend] [ghlgrv] [ghlphv] [SNRmin]"<<endl;
      return 0;
   }
 
/* ghlgrv & ghlphv
Gaussian halflength as percentage for group and phase dispersion shift
defines how the amplitude and SNRs are weighted
*/
   FILE *ff;
   char buff[500];
   int i, nsta, npth, nblk;
   struct STATION sta[NSTA];
   struct STAPAIR *spr = NULL;
/* read in sta.name, sta.lon, sta.lat from station list */
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
   if( nsta == 0 ) {
      cerr<<"Error(main): empty station list!"<<endl;
      exit(0);
   }

/* read in sta1, sta2, dispf, snrf pdisp_g pdisp_p from file list */
   if((ff=fopen(arg[2],"r"))==NULL) {
      cout<<"Cannot open file "<<arg[2]<<endl;
      return -1;
   }
   for(nblk=0,i=0;fgets(buff, 500, ff)!=NULL;) {
      if( nblk*BLKs <= i ) spr = (struct STAPAIR *) realloc (spr, (++nblk)*BLKs * sizeof(struct STAPAIR));
      if( sscanf(buff,"%s %s %s %s %s %s %f %s %s", spr[i].sta1, spr[i].sta2, spr[i].disp_pf, spr[i].snr_pf, spr[i].disp_nf, spr[i].snr_nf, &(spr[i].daynum), spr[i].pdisp_g, spr[i].pdisp_p) != 9 ) {
	 cerr<<"Warning(main): format error in file "<<arg[2]<<endl;
	 continue;
      }
      i++;
   }
   npth = i;
   fclose(ff);
   if( npth == 0 ) {
      cerr<<"Error(main): empty file list!"<<endl;
      exit(0);
   }

/* main loop. Process data from each path. Compute dist and snr, and as an indicator of FTAN quality,
   compute summation of SNR&amp weighted by distance between observed and predicted dispersion curves */
   int ipth, isp, isn, nlag;
   int flagp, flagn;
   float snrsig_pos, snrsig_neg, ampsig_pos, ampsig_neg;
   double dist, azi1, azi2;
   struct LAG_DATA datapos, dataneg;
   float TShiftG, sigmaG, TShiftP, sigmaP;
   perl = atof(arg[3]); perh = atof(arg[4]); 
   ghlgrv = atof(arg[5]); ghlphv = atof(arg[6]);
   snrmin = atof(arg[7]);
   alphagrv = -0.5/(ghlgrv*ghlgrv); alphaphv = -0.5/(ghlphv*ghlphv);

   sprintf(buff, "Disp_quality_%.1f_%.1f", perl, perh);
   ff = fopen(buff, "w");
   fprintf(ff, "sta1(1) lat1(2) lon1(3)  sta2(4) lat2(5) lon2(6)  dist(7) azi1(8) azi2(9) dnum(10) : snrsig_pos(12) ampsig_pos/dnum(13) snrsig_neg(14) ampsig_neg/dnum(15) GrvTimeShift(16) GrvUncertainty(17) PhvTimeShift(18) PhvUnvertainty(19)\n");
   //omp_init_lock(&readlock);
   //#pragma omp parallel for lastprivate(ipth) schedule (static,1)
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
   // compute weights and averaged SNR&amp for each lag
      nlag=0;
      if( AmpSNRSummation( spr[ipth].disp_pf, spr[ipth].snr_pf, spr[ipth].pdisp_g, spr[ipth].pdisp_p, &snrsig_pos, &ampsig_pos, &datapos ) ) nlag++;
      if( AmpSNRSummation( spr[ipth].disp_nf, spr[ipth].snr_nf, spr[ipth].pdisp_g, spr[ipth].pdisp_p, &snrsig_neg, &ampsig_neg, &dataneg ) ) nlag++;
      if( nlag == 0 ) {
	 cout<<"   Bad measurements in both lags. Skiped!"<<endl;
	 continue;
      }
   // compute averaged timeshift and its uncertainty
      if( nlag == 2 ) TimeShift( &datapos, &dataneg, dist, &TShiftG, &sigmaG, &TShiftP, &sigmaP );
      else { TShiftG = -1; sigmaG = -1; TShiftP = -1; sigmaP = -1; }
      cout<<"   Done"<<endl;
      //#pragma omp critical
      fprintf(ff, "%s %f %f  %s %f %f  %lf %lf %lf %f : %f %g  %f %g %f %f %f %f\n", sta[isp].name, sta[isp].lat, sta[isp].lon, sta[isn].name, sta[isn].lat, sta[isn].lon, dist, azi1, azi2, spr[ipth].daynum, snrsig_pos, ampsig_pos/spr[ipth].daynum, snrsig_neg, ampsig_neg/spr[ipth].daynum, TShiftG, sigmaG, TShiftP, sigmaP);
   }
   fclose(ff);
   //omp_destroy_lock(&readlock);

   free(spr);
 
   return 1;
}
