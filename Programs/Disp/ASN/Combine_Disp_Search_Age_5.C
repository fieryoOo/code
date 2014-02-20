#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include "/home/tianye/code/MyLib/Dis_Azi.h"
using namespace std;

#define NSTA 1000
#define NPTH 1000

//void calc_age(float lon1, float lat1, float lon2, float lat2, float **age, float lonmn, float latmn, int nlon, int nlat, float dgr, float *agemn, float *agemx, float *ageav){
void calc_age(float lon1, float lat1, float lon2, float lat2, char *agefile, float cc, float csw, float clw, float cse, float cle, float *ageav, float *std){
   char buff[300];
   float width=10;
   sprintf(buff, "project %s -C%f/%f -E%f/%f -Lw -Q -S -W%f/%f -Fzpxy > age.temp", agefile, lon1, lat1, lon2, lat2, -width, width);
   system(buff);
   FILE *ff = fopen("age.temp", "r");
   float age[5000], p[5000], lon[5000], lat[5000], dis[5000], vel, T;
   int i, n;
   double dist;
   calc_dist(lat1, lon1, lat2, lon2, &dist);
   for(i=0;;) {
      if(fgets(buff, 300, ff)==NULL)break;
      sscanf(buff,"%f %f %f %f",&age[i], &p[i], &lon[i], &lat[i]);
      if(age[i]==age[i]) i++;
   }
   fclose(ff);
   n=i;
   if(n==0) { *ageav = -1; return; }
   for(i=1;i<n-1;i++) dis[i] = (p[i+1]-p[i-1])/2.;
   dis[0] = (p[1]+p[0])/2.; dis[n-1] = dist-(p[n-1]+p[n-2])/2.;
   T=0;
   for(i=0;i<n;i++) {
      if(lat[i]>2.49*lon[i]-527.36) vel = cc+csw*sqrt(age[i])+clw*age[i];
      else vel = cc+cse*sqrt(age[i])+cle*age[i];
      T += dis[i]/vel;
   }
   vel = dist/T;
   float B = (2*cc/cle - 2*vel/cle - cse*cse/cle/cle), C = pow((vel-cc)/cle,2), D = B*B-4*C;
   if(D<0) {cout<<"negative!!"<<endl; exit(0);}
   float ageav1 = -B/2.+sqrt(D)/2., ageav2 = -B/2.-sqrt(D)/2.;
   float std1 = 0, std2 = 0;
   for(i=0;i<n;i++) { std1 = std1 + pow((age[i]-ageav1),2); std2 = std2 + pow((age[i]-ageav2),2); }
cout<<"1: "<<ageav1<<" "<<std1<<"  2: "<<ageav2<<" "<<std2<<endl;
   if(std1 < std2) { *ageav = ageav1; *std = sqrt(std1/(n-1)); }
   else { *ageav = ageav2; *std = sqrt(std2/(n-1)); }
}

float calc_diff(float phv, float lon1, float lat1, float lon2, float lat2, float cct, float cswt, float clwt, float cset, float clet, char *fname) {
   char buff[300];
   float width=10.;
   sprintf(buff, "project %s -C%f/%f -E%f/%f -Lw -Q -S -W%f/%f -Fzpxy > temp1", fname, lon1, lat1, lon2, lat2, -width, width);
   system(buff);
   FILE *ff = fopen("temp1", "r");
   float age[5000], p[5000], lon[5000], lat[5000], dis[5000], vel, T;
   int i, n;
   double dist;
   calc_dist(lat1, lon1, lat2, lon2, &dist);
   for(i=0;;) {
      if(fgets(buff, 300, ff)==NULL)break;
      sscanf(buff,"%f %f %f %f",&age[i], &p[i], &lon[i], &lat[i]);
      if(age[i]==age[i]) i++;
   }
   fclose(ff);
   n=i;
   if(n==0) return -1.;
   for(i=1;i<n-1;i++) dis[i] = (p[i+1]-p[i-1])/2.;
   dis[0] = (p[1]+p[0])/2.; dis[n-1] = dist-(p[n-1]+p[n-2])/2.;
   T=0;
   for(i=0;i<n;i++) {
      if(lat[i]>2.49*lon[i]-527.36) vel = cct+cswt*sqrt(age[i])+clwt*age[i];
      else vel = cct+cset*sqrt(age[i])+clet*age[i];
      T += dis[i]/vel;
   }
   vel = dist/T;
   return fabs(vel-phv);
}

void Search_ph_age(float *phv, float *lon1, float *lat1, float *lon2, float *lat2, int iper, float *cc, float *csw,  float *clw, float *cse, float *cle, char *fname) {
   float ccl=3.36, cch=3.37, cswl=0.25, cswh=0.252, clwl=-0.08, clwh=-0.09, csel=0.25, cseh=0.252, clel=-0.08, cleh=-0.09;
   float cct, cswt, clwt, cset, clet;
   int i, ne, icc, isw, ilw, ise, ile, nstep = 2;
   float rsd[nstep][nstep][nstep][nstep][nstep], rsdmin=100., rsdt;
   for(icc=0;icc<nstep;icc++) {
cout<<endl<<"icc: "<<icc<<endl;
      cct = ccl + (cch-ccl)/(nstep-1)*icc;
      for(isw=0;isw<nstep;isw++) {
cout<<endl<<"isw: "<<isw<<endl;
         cswt = cswl + (cswh-cswl)/(nstep-1)*isw;
         for(ilw=0;ilw<nstep;ilw++) {
            clwt = clwl + (clwh-clwl)/(nstep-1)*ilw;
            for(ise=0;ise<nstep;ise++) {
               cset = csel + (cseh-csel)/(nstep-1)*ise;
               for(ile=0;ile<nstep;ile++) {
                  clet = clel + (cleh-clel)/(nstep-1)*ile;
                  rsd[icc][isw][ilw][ise][ile] = 0;
                  for(ne=0,i=0;i<iper;i++) {
                     rsdt = calc_diff(phv[i], lon1[i], lat1[i], lon2[i], lat2[i], cct, cswt, clwt, cset, clet, fname);
                     if(rsdt<0) continue;
                     rsd[icc][isw][ilw][ise][ile] += rsdt*rsdt;
                     ne++;
                  }
                  rsd[icc][isw][ilw][ise][ile] = sqrt(rsd[icc][isw][ilw][ise][ile]/(ne-1));
            cout<<rsd[icc][isw][ilw][ise][ile]<<endl;
                  if(rsdmin>rsd[icc][isw][ilw][ise][ile]) {
                     rsdmin = rsd[icc][isw][ilw][ise][ile];
                     *cc = cct;
                     *csw = cswt;
                     *clw = clwt;
                     *cse = cset;
                     *cle = clet;
                  }
               }
            }
         }
      }
   }
}

void InserSort1(float *arr, float *dat, int n) {
   int i,j;
   float ftmp0, ftmp1;
   for(i=1;i<n;i++) {
      ftmp0=arr[i]; ftmp1=dat[i];
      for(j=i;j>0 && ftmp0<arr[j-1]; j--) { arr[j]=arr[j-1]; dat[j]=dat[j-1]; }
      arr[j]=ftmp0; dat[i]=ftmp1;
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

void get_snr(char *fname, float pero, float *snro) {
   FILE *ff;
   char buff[300];
   int i, n;
   float per[100], snr[100], ftmp;
   *snro = 0;
   if((ff=fopen(fname,"r"))==NULL) return;
   for(i=0;fgets(buff, 300, ff)!=NULL;i++)
      if(sscanf(buff, "%f %f %f", &per[i], &ftmp, &snr[i])!=3) continue;
   fclose(ff);
   n = i;
   InserSort1(per,snr,n);
   for(i=0;i<n && per[i]<pero; i++) {}
   if(i==0 || i==n) return;
   *snro = snr[i-1]+(pero-per[i-1])/(per[i]-per[i-1])*(snr[i]-snr[i-1]);
}

void get_vel(char *fname, float pero, float *grvo, float *phvo) {
   FILE *ff;
   char buff[300];
   int i, n, itmp;
   float per[100], grv[100], phv[100], ftmp;
   *grvo = -1; *phvo = -1;
   if((ff=fopen(fname,"r"))==NULL) return;
   for(i=0;fgets(buff, 300, ff)!=NULL;i++)
      if(sscanf(buff, "%d %f %f %f %f", &itmp, &ftmp, &per[i], &grv[i], &phv[i])!=5) continue;
   fclose(ff);
   n = i;
   InserSort2(per,grv,phv,n);
   for(i=0;i<n && per[i]<pero; i++) {}
   if(i==0 || i==n) return;
   *grvo = grv[i-1]+(pero-per[i-1])/(per[i]-per[i-1])*(grv[i]-grv[i-1]);
   *phvo = phv[i-1]+(pero-per[i-1])/(per[i]-per[i-1])*(phv[i]-phv[i-1]);
}

int ph_match(char *nfname, char *pfname, float pero) {
   int i;
   float per[3];
   float grvn[3], phvn[3], grvp[3], phvp[3];
   per[0] = pero/1.05, per[1] = pero, per[2] = pero/0.95;
   for(i=0;i<3;i++) {
      get_vel(nfname, per[i], &grvn[i], &phvn[i]);
      get_vel(pfname, per[i], &grvp[i], &phvp[i]);
      if(fabs(phvn[i]-phvp[i])>0.05) return 0;
   }
   return 1;
}

main(int na, char *arg[])
{
  if (na!=5)
    {
      cout<<"usage: "<<arg[0]<<" [station.lst] [sta1-sta2-Dispposf-SNRp-Dispnegf-SNRn-Dispsymf list] [Age_map (res 0.1)] [per]"<<endl;
      return 0;
    }
   
   FILE *ff, *fsnrp, *fsnrn, *fdspp, *fdspn;
   char buff[300], ctmp[100];
   char stan[NSTA][6], sta1[NPTH][6], sta2[NPTH][6];
   char dispsf[NPTH][100], disppf[NPTH][100], dispnf[NPTH][100], snrsf[NPTH][100], snrpf[NPTH][100], snrnf[NPTH][100];
   int i, j, iper, ipth, npth, nsta, itmp, flag;
   float slon[NSTA], slat[NSTA], ftmp, snrpos, snrneg;
   float lon1[NPTH], lon2[NPTH], lat1[NPTH], lat2[NPTH];
   float grv[NPTH], phv[NPTH], grvtmp, phvtmp, per[NPTH], pertmp;
   int ista1[NPTH], ista2[NPTH];
   float agemn, agemx, ageav, stdtmp, agemin[NPTH], agemax[NPTH], ageavg[NPTH], std[NPTH], snr[NPTH], distfact[NPTH], azii[NPTH];
   double dist, azi;
   float  pero=atof(arg[4]);


//read in stan, slon, slat from station list
   if((ff=fopen(arg[1],"r"))==NULL) {
      cout<<"Cannot open file "<<arg[1]<<endl;
      return -1;
   }
   for(i=0;;i++) {
      if(fgets(buff, 300, ff)==NULL) break;
      sscanf(buff,"%s %f %f", stan[i], &slon[i], &slat[i]);
      if(slon[i]<0) slon[i] += 360.;
   }
   nsta=i;
   fclose(ff);
//read in sta1, sta2, dispf, snrf from file list
   if((ff=fopen(arg[2],"r"))==NULL) {
      cout<<"Cannot open file "<<arg[2]<<endl;
      return -1;
   }
   for(i=0;;i++) {
      if(fgets(buff, 300, ff)==NULL) break;
      sscanf(buff,"%s %s %s %s %s %s %s", sta1[i], sta2[i], disppf[i], snrpf[i], dispnf[i], snrnf[i], dispsf[i]);
   }
   npth = i;
   fclose(ff);
//read in dispersion data from each path, check for dist and snr
   for(iper=0,ipth=0;ipth<npth;ipth++) {
      for(i=0;i<nsta;i++) if(strcmp(sta1[ipth],stan[i])==0) break;
      if(i==nsta) continue;
      lon1[iper] = slon[i]; lat1[iper] = slat[i]; ista1[iper] = i;
      for(i=0;i<nsta;i++) if(strcmp(sta2[ipth],stan[i])==0) break;
      if(i==nsta) continue;
      lon2[iper] = slon[i]; lat2[iper] = slat[i]; ista2[iper] = i;
      //discard path according to locations
      if(lat1[iper]>2.49*lon1[iper]-528 && lat2[iper]>2.49*lon2[iper]-528) continue;
      calc_dist(lat1[iper], lon1[iper], lat2[iper], lon2[iper], &dist);
      distfact[iper] = dist/pero;
      if(distfact[iper] < 12.) continue;
      calc_azimuth(lat1[iper], lon1[iper], lat2[iper], lon2[iper], &azi);
      azii[iper] = azi;
      cout<<"Working on path "<<sta1[ipth]<<" - "<<sta2[ipth]<<endl;
      get_snr(snrpf[ipth], pero, &snrpos);
      get_snr(snrnf[ipth], pero, &snrneg);
      if(snrpos < 5 && snrneg < 5) continue;
      flag = 0;
      if(snrpos > 5 && snrneg > 5 && ph_match(dispnf[ipth], disppf[ipth], pero)) {
            snr[iper] = (snrpos+snrneg)/1.41;
            get_vel(dispsf[ipth], pero, &grv[iper], &phv[iper]);
            if(grv[iper]<0) continue;
            flag = 1;
      }
      else if(snrpos > snrneg) { get_vel(disppf[ipth], pero, &grv[iper], &phv[iper]); snr[iper] = snrpos; }
      else { get_vel(dispnf[ipth], pero, &grv[iper], &phv[iper]); snr[iper] = snrneg; }
      //calc_age(lon1, lat1, lon2, lat2, arg[3], &agemin[iper], &agemax[iper], &ageavg[iper], &std[iper]);
      iper++;
   }
cout<<"npath: "<<iper<<endl;
  
   float cc, csw, clw, cse, cle;
   Search_ph_age(phv, lon1, lat1, lon2, lat2, iper, &cc, &csw, &clw, &cse, &cle, arg[3]);
   for(i=0;i<iper;i++) calc_age(lon1[i], lat1[i], lon2[i], lat2[i], arg[3], cc, csw, clw, cse, cle, &ageavg[i], &std[i]);
   cout<<cc<<" "<<csw<<" "<<clw<<" "<<cse<<" "<<cle<<endl;

/*
   sprintf(buff, "Age_Grv_%.1fsec",pero);
   ff=fopen(buff, "w");
   for(i=0;i<iper;i++) fprintf(ff, "%f %f %f %f %f %f %s %f %f %s %f %f %f %f %f\n", pero, grv[i], ageavg[i], std[i], snr[i], distfact[i], stan[ista1[i]], slon[ista1[i]], slat[ista1[i]], stan[ista2[i]], slon[ista2[i]], slat[ista2[i]], azii[i], agemin[i], agemax[i]);
   fclose(ff);
*/
   sprintf(buff, "Age_Phv_%.1fsec_GS", pero);
   ff=fopen(buff, "w");
   for(i=0;i<iper;i++) fprintf(ff, "%f %f %f %f %f %f %s %f %f %s %f %f %f\n", pero, phv[i], ageavg[i], std[i], snr[i], distfact[i], stan[ista1[i]], slon[ista1[i]], slat[ista1[i]], stan[ista2[i]], slon[ista2[i]], slat[ista2[i]], azii[i]);
   fclose(ff);

  return 1;
}
