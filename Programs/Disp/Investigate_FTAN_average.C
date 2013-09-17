/*
This code computes the averaged dispersion information around the strongest period in a given period band.
No quality control other than the SNR critera is performed. Data points with a group velocity out of 1.5 
standard deviation are discarded
Output: (distance, azimuth, averaged snr, averaged amp/dnum, averaged grv, and averaged phv of both positive 
and negative lags)
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
#define BLKs 5000
#define NFRQ 100

float perl, perh, snrmin;

struct STATION {
   char name[6];
   float lon, lat;
};

struct STA_PAIR {
   char sta1[6], sta2[6];
   char disp_pf[150], disp_nf[150];
   char snr_pf[150], snr_nf[150];
   float daynum;
};

struct SPR_DATA {
   int npt;
   float per[NFRQ], snr[NFRQ], wsum;
   float grv[NFRQ], phv[NFRQ], amp[NFRQ];
};

void InserSort(struct SPR_DATA *dat) {
   int i,j;
   float ftmp[5];
   for(i=1;i<dat->npt;i++) {
      ftmp[0]=dat->per[i]; ftmp[1]=dat->snr[i];
      ftmp[2]=dat->grv[i]; ftmp[3]=dat->phv[i]; ftmp[4]=dat->amp[i];
      for(j=i;j>0 && ftmp[0]<dat->per[j-1]; j--) { 
	 dat->per[j]=dat->per[j-1]; dat->snr[j]=dat->snr[j-1]; 
	 dat->grv[j]=dat->grv[j-1]; dat->phv[j]=dat->phv[j-1];
	 dat->amp[j]=dat->amp[j-1];
      }
      dat->per[j]=ftmp[0]; dat->snr[j]=ftmp[1];
      dat->grv[j]=ftmp[2]; dat->phv[j]=ftmp[3];
      dat->amp[j]=ftmp[4];
   }
}

float Avg( float *arr, int npt ) {
   float sum=0.;
   for(int i=0; i<npt; i++) sum += arr[i];
   return sum/npt;
}
int DataSelect( char *fdisp, char *famp, struct SPR_DATA *sdata ) {
   FILE *fin;
   char buff[300];
   int i, ief, itmp;
   float ftmp;
   struct SPR_DATA datatmp;
   // read fdisp
   if((fin=fopen(fdisp,"r"))==NULL) return 0;
   for(i=0;fgets(buff, 300, fin)!=NULL;) {
      if(sscanf( buff, "%d %f %f %f %f", &itmp, &ftmp, &(datatmp.per[i]), &(datatmp.grv[i]), &(datatmp.phv[i]) )!=5) continue;
      if( datatmp.per[i]<perl || datatmp.per[i]>perh ) continue;
      i++;
   }
   fclose(fin);
   // read famp
   if((fin=fopen(famp,"r"))==NULL) return 0;
   for(i=0;fgets(buff, 300, fin)!=NULL;) {
      if(sscanf( buff, "%f %f %f", &ftmp, &(datatmp.amp[i]), &(datatmp.snr[i]) )!=3) continue;
      if( ftmp<perl || ftmp>perh ) continue;
      if( ftmp != datatmp.per[i] ) {
	 cerr<<"Warning: mismatch between file "<<fdisp<<" and "<<famp<<endl;
	 return 0;
      }
      i++;
   }
   itmp = i;
   fclose(fin);
   // sort by period
   InserSort(&datatmp);
   // compute avg and std of vel
   float vstd, vavg = Avg(datatmp.grv, itmp);
   for(vstd=0.,i=0;i<itmp;i++) {
      ftmp = datatmp.grv[i] - vavg;
      vstd += ftmp*ftmp;
   }
   vstd = 1.5*sqrt(vstd/(itmp-1));
   // discard by vstd/snrmin and assign weights
   for(ief=0,i=0;i<itmp;i++) {
      if( fabs(datatmp.grv[i]-vavg) > vstd ) continue;
      if( datatmp.snr[i] < snrmin ) continue;
      if( i > ief ) {
	 ftmp = datatmp.snr[i] * datatmp.snr[i];
	 sdata->per[ief] = datatmp.per[i]*ftmp; sdata->snr[ief] = datatmp.snr[i]*ftmp;
	 sdata->grv[ief] = datatmp.grv[i]*ftmp; sdata->phv[ief] = datatmp.phv[i]*ftmp;
	 sdata->amp[ief] = datatmp.amp[i]*ftmp; sdata->wsum += ftmp;
      }
      ief++;
   }
   sdata->npt = ief;
   return 1;  
}


main(int na, char *arg[])
{
   if(na!=6) {
      cout<<"usage: "<<arg[0]<<" [station.lst] [sta1-sta2-Dispposf-SNRp-Dispnegf-SNRn-daynum list] [period_lowerend] [period_higherend] [min snr]"<<endl;
      return 0;
   }
 
   FILE *ff;
   char buff[300];
   int i, nsta, npth, nblk;
   struct STATION sta[NSTA];
   struct STA_PAIR *spr = NULL;
// read in sta.name, sta.lon, sta.lat from station list
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

// read in sta1, sta2, dispf, snrf from file list
   if((ff=fopen(arg[2],"r"))==NULL) {
      cout<<"Cannot open file "<<arg[2]<<endl;
      return -1;
   }
   for(nblk=0,i=0;fgets(buff, 300, ff)!=NULL;i++) {
      if( nblk*BLKs <= i ) spr = (struct STA_PAIR *) realloc (spr, (++nblk)*BLKs * sizeof(struct STA_PAIR));
      sscanf(buff,"%s %s %s %s %s %s %f", spr[i].sta1, spr[i].sta2, spr[i].disp_pf, spr[i].snr_pf, spr[i].disp_nf, spr[i].snr_nf, &(spr[i].daynum));
   }
   npth = i;
   fclose(ff);
// read in parameters
   perl = atof(arg[3]); perh = atof(arg[4]);
   snrmin = atof(arg[5]);
// read in dispersion data from each path, check for dist and snr
   int ipth, isp, isn, np, nn;
   int flagp, flagn;
   double dist, azi1, azi2;
   float perp, pern, grvp, grvn, phvp, phvn, ampp, ampn, snrp, snrn;
   struct SPR_DATA datapos, dataneg;
   sprintf(buff, "Disp_info_avg_%.1fto%.1fsec", perl, perh);
   ff = fopen(buff, "w");
   fprintf(ff, "sta1 lat1 lon1  sta2 lat2 lon2  dist azi1 azi2 dnum : per-appa snrp ampp/dnum grvp phvp  per-appa snrn ampn/dnum grvn phvn\n");
   for(ipth=0;ipth<npth;ipth++) {
      // search for sta1 and sta2 in station list
      for(isp=0;isp<nsta;isp++) if(strcmp(spr[ipth].sta1, sta[isp].name)==0) break;
      if(isp==nsta) continue;
      for(isn=0;isn<nsta;isn++) if(strcmp(spr[ipth].sta2, sta[isn].name)==0) break;
      if(isn==nsta) continue;
      cout<<"Extracting information for path "<<spr[ipth].sta1<<" - "<<spr[ipth].sta2<<endl;
      calc_dist(sta[isp].lat, sta[isp].lon, sta[isn].lat, sta[isn].lon, &dist);
      calc_azimuth(sta[isp].lat, sta[isp].lon, sta[isn].lat, sta[isn].lon, &azi1);
      calc_azimuth(sta[isn].lat, sta[isn].lon, sta[isp].lat, sta[isp].lon, &azi2);
      // read in the disp amd amp_snr file, discard data points by grv-std and snrmin and weight data by energy (snr sqaure)
      if( ! DataSelect( spr[ipth].disp_pf, spr[ipth].snr_pf, &datapos ) ) continue;
      if( ! DataSelect( spr[ipth].disp_nf, spr[ipth].snr_nf, &dataneg ) ) continue;
      // get averaged per, amp, snr, grv, and phv
      np = datapos.npt; nn = dataneg.npt; //npt
      perp = Avg(datapos.per, np); pern = Avg(dataneg.per, nn); //per
      snrp = Avg(datapos.snr, np); snrn = Avg(dataneg.snr, nn); //snr
      ampp = Avg(datapos.amp, np); ampn = Avg(dataneg.amp, nn); //amp
      grvp = Avg(datapos.grv, np); grvn = Avg(dataneg.grv, nn); //grv
      phvp = Avg(datapos.phv, np); phvn = Avg(dataneg.phv, nn); //phv
      // output
      fprintf(ff, "%s %f %f  %s %f %f  %lf %lf %lf %f : %f %f %g %f %f  %f %f %g %f %f\n", sta[isp].name, sta[isp].lat, sta[isp].lon, sta[isn].name, sta[isn].lat, sta[isn].lon, dist, azi1, azi2, spr[ipth].daynum, perp, snrp, ampp/spr[ipth].daynum, grvp, phvp, pern, snrn, ampn/spr[ipth].daynum, grvn, phvn);
   }
   fclose(ff);

   free(spr);
  
   return 1;
}
