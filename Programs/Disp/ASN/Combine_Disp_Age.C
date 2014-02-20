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
void calc_age(float lon1, float lat1, float lon2, float lat2, char *agefile, float *agemn, float *agemx, float *ageav, float *std){
   char buff[300];
   float width=10;
   sprintf(buff, "project %s -C%f/%f -E%f/%f -Lw -Q -S -W%f/%f -Fzp > age.temp", agefile, lon1, lat1, lon2, lat2, -width, width);
   system(buff);
   FILE *ff = fopen("age.temp", "r");
   float age[5000], p[5000], dis;
   int i, n;
   *agemn=99; *agemx=0; *ageav=0; *std=0;
   for(i=0;;i++) {
      if(fgets(buff, 300, ff)==NULL)break;
      sscanf(buff,"%f %f",&age[i], &p[i]);
      if(*agemn>age[i]) *agemn = age[i];
      if(*agemx<age[i]) *agemx = age[i];
   }
   fclose(ff);
   n=i; dis=0.;
   for(i=1;i<n-1;i++) {
      //*ageav = *ageav+(p[i+1]-p[i-1])/2./(-0.01*pow(age[i]-2.5,2)+0.05*age[i]+3.68);
      *ageav = *ageav+(p[i+1]-p[i-1])/2.*age[i];
      dis += (p[i+1]-p[i-1])/2.;
   }
   //*ageav = dis/(*ageav); if(*ageav>3.867) *ageav = 3.867; *ageav = 5.-5.*sqrt(15.47-4*(*ageav));
   *ageav = (*ageav)/dis;
   for(i=0;i<n;i++) {
      *std = *std + pow((age[i]-*ageav),2);
   }
   *std = sqrt(*std/(n-1));
}

main(int na, char *arg[])
{
  if (na!=5)
    {
      cout<<"usage: "<<arg[0]<<" [station.lst] [sta1-sta2-Dispf-SNRf list] [Age_map (res 0.1)] [out_name]"<<endl;
      return 0;
    }
   
   FILE *ff, *fsnr, *fdsp;
   char buff[300], ctmp[100];
   char stan[NSTA][6], sta1[NPTH][6], sta2[NPTH][6], dispf[NPTH][100], snrf[NPTH][100];
   int i, j, iper, ipth, npth, nsta, itmp;
   float slon[NSTA], slat[NSTA], ftmp, snrpos, snrneg;
   float lon1, lon2, lat1, lat2;
   float grv[NPTH*20], phv[NPTH*20], per[NPTH*20], pertmp;
   int ista1[NPTH*20], ista2[NPTH*20];
   float agemn, agemx, ageav, stdtmp, agemin[NPTH*20], agemax[NPTH*20], ageavg[NPTH*20], std[NPTH*20], snr[NPTH*20], distfact[NPTH*20], azii[NPTH*20];
   double dist, azi;

//read in Age map
/*
   float lonmx=0., lonmn=360., latmn=90., latmx=-90;
   float lontmp, lattmp;
   if((ff=fopen(arg[3],"r"))==NULL) {
      cout<<"Cannot open file "<<arg[3]<<endl;
      return -1;
   }
   for(i=0;;i++) {
      if(fgets(buff, 300, ff)==NULL) break;
      sscanf(buff,"%f %f %f", &lontmp, &lattmp, &ftmp);
      if(lontmp<0) lontmp += 360.;
      if(lonmx<lontmp) lonmx=lontmp;
      if(lonmn>lontmp) lonmn=lontmp;
      if(latmx<lattmp) latmx=lattmp;
      if(latmn>lattmp) latmn=lattmp;
   }
   rewind(ff);
   float dgr=0.1;
   int nlon = int((lonmx-lonmn)/dgr+0.5)+1, nlat = int((latmx-latmn)/dgr+0.5)+1;
   int ilon, ilat;
   float **age;
   age = (float **) malloc (nlon * sizeof(double *));
   for(i=0;i<nlon;i++) age[i] = (float *) malloc (nlat * sizeof(double *));
   float agemn, agemx, ageav, agemin[NPTH*20], agemax[NPTH*20], ageavg[NPTH*20];
   for(i=0;;i++) {
      if(fgets(buff, 300, ff)==NULL) break;
      if(sscanf(buff,"%f %f %f", &lontmp, &lattmp, &ftmp)!=3) continue;
      if(lontmp<0) lontmp += 360.;
      ilon = int((lontmp-lonmn)/dgr+0.5);
      ilat = int((lattmp-latmn)/dgr+0.5);
      age[ilon][ilat] = ftmp;
   }
   fclose(ff);
*/ 

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
      sscanf(buff,"%s %s %s %s", sta1[i], sta2[i], dispf[i], snrf[i]);
   }
   npth = i;
   fclose(ff);
//read in dispersion data from each path, check for dist and snr
   for(iper=0,ipth=0;ipth<npth;ipth++) {
      for(i=0;i<nsta;i++) if(strcmp(sta1[ipth],stan[i])==0) break;
      if(i==nsta) continue;
      lon1 = slon[i]; lat1 = slat[i]; ista1[iper] = i;
      for(i=0;i<nsta;i++) if(strcmp(sta2[ipth],stan[i])==0) break;
      if(i==nsta) continue;
      lon2 = slon[i]; lat2 = slat[i]; ista2[iper] = i;
      calc_dist(lat1, lon1, lat2, lon2, &dist);
      calc_azimuth(lat1, lon1, lat2, lon2, &azi);
      if((fsnr=fopen(snrf[ipth],"r"))==NULL) {
         cout<<"Cannot open file "<<snrf[ipth]<<endl;
         continue;
      }
      if((fdsp=fopen(dispf[ipth],"r"))==NULL) {
         cout<<"Cannot open file "<<dispf[ipth]<<endl;
         continue;
      }
      cout<<"Working on file "<<dispf[ipth]<<endl;
      //calc_age(lon1, lat1, lon2, lat2, age, lonmn, latmn, nlon, nlat, dgr, &agemn, &agemx, &ageav);
      calc_age(lon1, lat1, lon2, lat2, arg[3], &agemn, &agemx, &ageav, &stdtmp);
      for(i=0;;i++) {
         if(fgets(buff, 300, fsnr)==NULL) break;
         if(sscanf(buff,"%f %f %f %f %f", &pertmp, &ftmp, &snrpos, &ftmp, &snrneg)!=5) {cout<<"Wrong format in file "<<snrf[ipth]<<endl; return -1;}
         if(fgets(buff, 300, fdsp)==NULL) break;
         if(sscanf(buff,"%d %f %f %f %f", &itmp, &ftmp, &per[iper], &grv[iper], &phv[iper])!=5) {cout<<"Wrong format in file "<<dispf[ipth]<<endl; return -1;}
         if(pertmp!=per[iper]) {cout<<"Mismatch between disp and snr file"<<endl; break;}
         if(per[iper]*10. > dist) continue;
         if( snrpos+snrneg < 12 ) continue;
         agemin[iper] = agemn; agemax[iper] = agemx; ageavg[iper] = ageav; std[iper] = stdtmp;
         snr[iper] = snrpos+snrneg; distfact[iper] = dist/per[iper]; azii[iper] = azi;
         iper++;
         ista1[iper] = ista1[iper-1];
         ista2[iper] = ista2[iper-1];
      }
      fclose(fdsp); fclose(fsnr);
   }
cout<<"nper: "<<iper<<endl;
  
   sprintf(buff, "%s.gr_age", arg[4]);
   ff=fopen(buff, "w");
   for(i=0;i<iper;i++) fprintf(ff, "%f %f %f %f %f %f %s %f %f %s %f %f %f %f %f\n", per[i], grv[i], ageavg[i], std[i], snr[i], distfact[i], stan[ista1[i]], slon[ista1[i]], slat[ista1[i]], stan[ista2[i]], slon[ista2[i]], slat[ista2[i]], azii[i], agemin[i], agemax[i]);
   fclose(ff);
   sprintf(buff, "%s.ph_age", arg[4]);
   ff=fopen(buff, "w");
   for(i=0;i<iper;i++) fprintf(ff, "%f %f %f %f %f %f %s %f %f %s %f %f %f %f %f\n", per[i], phv[i], ageavg[i], std[i], snr[i], distfact[i], stan[ista1[i]], slon[ista1[i]], slat[ista1[i]], stan[ista2[i]], slon[ista2[i]], slat[ista2[i]], azii[i], agemin[i], agemax[i]);
   fclose(ff);

  return 1;
}
