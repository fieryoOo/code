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
#define NPER 100

//least square fit
void fit_line (float *datx, float *daty, int ndat, int indep, float *a, float *b, float *std) {
   FILE *ff;
   char buff[300];
   int i;
   float X=0, Y=0, X2=0, Y2=0, XY=0;
   for(i=0;i<ndat;i++) {
      X += datx[i];
      Y += daty[i];
      X2 += pow(datx[i],2);
      Y2 += pow(daty[i],2);
      XY += datx[i]*daty[i];
   }
   if(indep==0) {
      *a = (ndat*XY-X*Y)/(ndat*X2-X*X);
      *b = (-X*XY+X2*Y)/(ndat*X2-X*X);
   }
   else if(indep==1) {
      *a=(ndat*Y2-Y*Y)/(ndat*XY-X*Y);
      *b=(Y*XY-Y2*X)/(ndat*XY-X*Y);
   }
   else {
      cout<<"Line_fit: Wrong input for indep var, stopped!"<<endl;
      exit(0);
   }
   *std=0;
   for(i=0;i<ndat;i++) *std += pow(daty[i]-(*a)*datx[i]-(*b),2);
   if(indep==1) *std = *std/(*a)/(*a);
   *std = sqrt(*std/(ndat-1));
}

void fit_line_2(char *fname, int indep, float *a, float *b, float *std) {
   FILE *ff;
   char buff[300];
   int i, j, ndat;
   float datx[NPTH], daty[NPTH];
   if((ff=fopen(fname,"r"))==NULL) {
      cout<<"Cannot open file "<<fname<<endl;
      exit(0);
   }
   for(i=0;;i++) {
      if(fgets(buff,300,ff)==NULL) break;
      sscanf(buff,"%f %f",&datx[i],&daty[i]);
   }
   fclose(ff);
   ndat=i;
  
   fit_line( datx, daty, ndat, 0, a, b, std );
   for(i=0;;i++) {
      if(i>=ndat) break;
      if(fabs(daty[i]-datx[i]*(*a)-(*b))<(*std)*2) continue;
      for(j=i;j<ndat-1;j++) { datx[j] = datx[j+1]; daty[j] = daty[j+1]; }
      ndat--;
   }
   fit_line( datx, daty, ndat, 0, a, b, std );
}

//void calc_age(float lon1, float lat1, float lon2, float lat2, float **age, float lonmn, float latmn, int nlon, int nlat, float dgr, float *agemn, float *agemx, float *ageav){
void calc_age(float lon1, float lat1, float lon2, float lat2, float dist, char *agefile, float a[], float b[], int nper, float *agemn, float *agemx, float *ageav, float *std, char *outname){
   char buff[300], name[100];
   float width=10;
   sprintf(name, "%s.age_dis", outname);
   sprintf(buff, "project %s -C%f/%f -E%f/%f -Lw -Q -S -W%f/%f -Fzp > %s", agefile, lon1, lat1, lon2, lat2, -width, width, name);
   system(buff);
   FILE *ff = fopen(name, "r");
   float age[5000], p[5000], dis[5000], vel, T;
   int i, iper, n;
   *agemn=99; *agemx=-1;
   for(i=0;;) {
      if(fgets(buff, 300, ff)==NULL)break;
      sscanf(buff,"%f %f",&age[i], &p[i]);
      if(age[i]==age[i]) {
         if(*agemn>age[i]) *agemn = age[i];
         if(*agemx<age[i]) *agemx = age[i];
         i++;
      }
   }
   fclose(ff);
   n=i;
   if(n==0) return;
   for(i=1;i<n-1;i++) dis[i] = (p[i+1]-p[i-1])/2.;
   dis[0] = (p[1]+p[0])/2.; dis[n-1] = dist-(p[n-1]+p[n-2])/2.;
   for(iper=0; iper<nper; iper++) {
      ageav[iper]=0; *std=0; T=0;
      for(i=0;i<n;i++) {
         vel = (age[i] - b[iper])/a[iper];
         T += dis[i]/vel;
      }
      vel = dist/T;
      ageav[iper] = a[iper]*vel+b[iper];
 cout<<lon1<<" "<<lat1<<" "<<lon2<<" "<<lat2<<" "<<vel<<" "<<T<<" "<<ageav[iper]<<endl;
   foreach per (6.0 6.5 7.0 7.5 8.0 8.5 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 21.0 22.0 23.0 24.0 25.0)
      for(i=0;i<n;i++) *std = *std + pow((age[i]-ageav[iper]),2);
      *std = sqrt(*std/(n-1));
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


main(int na, char *arg[])
{
  if (na!=5)
    {
      cout<<"usage: "<<arg[0]<<" [station.lst] [sta1-sta2 list] [Age_map (res 0.1)] [Age_Phv.flst (per file)]"<<endl;
      return 0;
    }
   
   FILE *ff, *fsnrp, *fsnrn, *fdspp, *fdspn;
   char buff[300], name[100], ctmp[100];
   char stan[NSTA][6], sta1[NPTH][6], sta2[NPTH][6];
   int i, j, iper, nper, ipth, npth, npthe, nsta, itmp, flag;
   float slon[NSTA], slat[NSTA], ftmp;
   float lon1, lon2, lat1, lat2;
   float per[NPER], a[NPER], b[NPER], stdtmp;
   int ista1[NPTH], ista2[NPTH];
   float agemn, agemx, ageav, agemin[NPTH], agemax[NPTH], ageavg[NPTH][NPER], std[NPTH], distfact[NPTH], azii[NPTH];
   double dist, azi;

//read in Age_Phv file names at each period
   if((ff=fopen(arg[4], "r"))==NULL) {
      cout<<"Cannot open file "<<arg[4]<<endl;
      return -1;
   }
   for(iper=0;;iper++) {
      if(fgets(buff, 300, ff)==NULL) break;
      sscanf(buff, "%f %s", &per[iper], &name[0]);
      //compute a[] b[] at each per[]
      fit_line_2(name, 0, &a[iper], &b[iper], &stdtmp);
cout<<per[iper]<<"sec: "<<a[iper]<<" "<<b[iper]<<" std: "<<stdtmp<<endl;
   }
   nper = iper;
   fclose(ff);
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
//read in sta1, sta2, path list
   if((ff=fopen(arg[2],"r"))==NULL) {
      cout<<"Cannot open file "<<arg[2]<<endl;
      return -1;
   }
   for(i=0;;i++) {
      if(fgets(buff, 300, ff)==NULL) break;
      sscanf(buff,"%s %s", sta1[i], sta2[i]);
   }
   npth = i;
   fclose(ff);
//read in dispersion data from each path, check for dist and snr
   for(npthe=0,ipth=0;ipth<npth;ipth++) {
      for(i=0;i<nsta;i++) if(strcmp(sta1[ipth],stan[i])==0) break;
      if(i==nsta) continue;
      lon1 = slon[i]; lat1 = slat[i]; ista1[npthe] = i;
      for(i=0;i<nsta;i++) if(strcmp(sta2[ipth],stan[i])==0) break;
      if(i==nsta) continue;
      lon2 = slon[i]; lat2 = slat[i]; ista2[npthe] = i;
      calc_dist(lat1, lon1, lat2, lon2, &dist);
      distfact[npthe] = dist;
      calc_azimuth(lat1, lon1, lat2, lon2, &azi);
      azii[npthe] = azi;
      cout<<"Working on path "<<sta1[ipth]<<" - "<<sta2[ipth]<<endl;
      sprintf(name, "%s_%s", sta1[ipth], sta2[ipth]);
      calc_age(lon1, lat1, lon2, lat2, distfact[npthe], arg[3], a, b, nper, &agemin[npthe], &agemax[npthe], &(ageavg[npthe][0]), &std[npthe], name);
      npthe++;
   }
cout<<"# of path processed: "<<npthe<<endl;
  
   sprintf(buff, "Path_Age.txt");
   ff=fopen(buff, "w");
   for(i=0;i<npthe;i++) if(ageavg[i]>0) fprintf(ff, "%s %s %f %f %f %f %f\n", stan[ista1[i]], stan[ista2[i]], ageavg[i], std[i], agemin[i], agemax[i], distfact[i]);
   fclose(ff);

  return 1;
}
