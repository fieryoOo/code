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

void calc_age(float lon1, float lat1, float lon2, float lat2, char *agefile, float cc, float cs, float cl, float *ageav, float *std){
   char buff[300];
   float width=8.;
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
      vel = cc+cs*sqrt(age[i])+cl*age[i];
      T += dis[i]/vel;
   }
   vel = dist/T;
   float B = (2*cc/cl - 2*vel/cl - cs*cs/cl/cl), C = pow((vel-cc)/cl,2), D = B*B-4*C;
   if(D<0) {cout<<"negative!!"<<endl; exit(0);}
   float ageav1 = -B/2.+sqrt(D)/2., ageav2 = -B/2.-sqrt(D)/2.;
   float std1 = 0, std2 = 0;
   for(i=0;i<n;i++) { std1 = std1 + pow((age[i]-ageav1),2); std2 = std2 + pow((age[i]-ageav2),2); }
// cout<<"1: "<<ageav1<<" "<<std1<<"  2: "<<ageav2<<" "<<std2<<endl;
   if(std1 < std2) { *ageav = ageav1; *std = sqrt(std1/(n-1)); }
   else { *ageav = ageav2; *std = sqrt(std2/(n-1)); }
}

void get_path_age (float lon1, float lat1, float lon2, float lat2, float age[], float dis[], int *npt, float *dist, char *fname) {
   char buff[300];
   float width=8., p[1000];
   sprintf(buff, "project %s -C%f/%f -E%f/%f -Lw -Q -S -W%f/%f -Fzp > temp1", fname, lon1, lat1, lon2, lat2, -width, width);
   system(buff);
   FILE *ff = fopen("temp1", "r");
   int i;
   double distmp;
   calc_dist(lat1, lon1, lat2, lon2, &distmp);
   *dist = distmp;
   for(i=0;;) {
      if(fgets(buff, 300, ff)==NULL)break;
      sscanf(buff,"%f %f",&age[i], &p[i]);
      if(age[i]==age[i]) i++;
   }
   fclose(ff);
   *npt=i;
   if(i==0) return;
   for(i=1;i<*npt-1;i++) dis[i] = (p[i+1]-p[i-1])/2.;
   dis[0] = (p[1]+p[0])/2.; dis[*npt-1] = *dist-(p[*npt-1]+p[*npt-2])/2.;
}

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

void Search_gr_age(float per, float *grv, float *lon1, float *lat1, float *lon2, float *lat2, int iper, float *cc, float *cs,  float *cl, int *flag, char *fname) {
   cout<<"Grid Searching for group vel..."<<endl;
   float ccl, cch, csl, csh, cll, clh;
   if(per==6) { ccl=0.5, cch=1.5, csl=-.5, csh=.5, cll=-0.15, clh=0.05; }
   else if(per==6.5) { ccl=0.5, cch=1.5, csl=-0.7, csh=0.3, cll=-0.1, clh=0.1; }
   else if(per==7.0) { ccl=1., cch=2., csl=-1., csh=0., cll=0., clh=0.2; }
   else if(per==7.5) { ccl=1.5, cch=3.5, csl=-2., csh=-0., cll=-0.5, clh=0.5; }
   else if(per==8.0) { ccl=2.5, cch=3.5, csl=-2., csh=-1., cll=0.15, clh=0.35; }
   else if(per==8.5) { ccl=2.5, cch=4.5, csl=-2., csh=0., cll=-0.5, clh=0.5; }
   else if(per==9.0) { ccl=2.5, cch=4.5, csl=-1.5, csh=0.5, cll=-0.5, clh=0.5; }
   else if(per==10.0) { ccl=2.8, cch=3.8, csl=-0.7, csh=0.3, cll=-0.1, clh=0.1; }
   else if(per==11.0) { ccl=2.5, cch=4.5, csl=-1., csh=1., cll=-0.5, clh=0.5; }
//   else if(per==12.0) { ccl=2.8, cch=3.8, csl=-0.2, csh=0.8, cll=-0.2, clh=0.; }
   else if(per==12.0) { ccl=2.6, cch=3.6, csl=0., csh=1., cll=-0.2, clh=-0.1; }
   else if(per==13.0) { ccl=3., cch=4., csl=-0.7, csh=0.3, cll=0.1, clh=0.3; }
   else if(per==14.0) { ccl=2.8, cch=3.8, csl=-.3, csh=.7, cll=-0.1, clh=0.1; }
   else if(per==15.0) { ccl=2.7, cch=3.7, csl=0., csh=1., cll=-0.2, clh=0.; }
   else if(per==16.0) { ccl=3.3, cch=4.3, csl=-1., csh=0., cll=0.18, clh=0.28; }//
   else if(per==17.0) { ccl=2., cch=4., csl=-0.5, csh=1.5, cll=-0.5, clh=0.5; }
   else if(per==18.0) { ccl=3.3, cch=4.3, csl=-1., csh=0., cll=0.15, clh=0.35; }
   else if(per==19.0) { ccl=2., cch=4., csl=-0.5, csh=1.5, cll=-0.5, clh=0.5; }
   else if(per==20.0) { ccl=1., cch=3., csl=1., csh=3., cll=-1.5, clh=0.5; }
   else { ccl=2., cch=4., csl=-1., csh=1., cll=-0.5, clh=0.5; }

   float cct, cst, clt;
   int i, ipt, ne, icc, ics, icl, nstep = 100, npt[iper];
   float *age[iper], *dis[iper], dist[iper];
   for(i=0;i<iper;i++) {age[i] = (float *) malloc (1000 * sizeof(float)); dis[i] = (float *) malloc (1000 * sizeof(float));}
//extract dist and age info for each path
   for(i=0;i<iper;i++) get_path_age(lon1[i], lat1[i], lon2[i], lat2[i], age[i], dis[i], &npt[i], &dist[i], fname);

   float rsd1[nstep][nstep][nstep][iper], rsd2[nstep][nstep][nstep][iper], rsdmin=100., rsdt, rsdm, rstd, T, vel;
   char name[100];
   sprintf(name,"rsdg_%.1f.txt",per);
   FILE *ff = fopen(name,"w");
   for(icc=0;icc<nstep;icc++) {
cout<<icc<<"/"<<nstep<<endl;
      cct = ccl + (cch-ccl)/(nstep-1)*icc;
      for(ics=0;ics<nstep;ics++) {
         cst = csl + (csh-csl)/(nstep-1)*ics;
         for(icl=0;icl<nstep;icl++) {
            clt = cll + (clh-cll)/(nstep-1)*icl;
          //calc rsd
            for(i=0;i<iper;i++) {
               if(npt[i]==0) { rsd1[icc][ics][icl][i] = -1.; continue; }
               T=0.;
               for(ipt=0;ipt<npt[i];ipt++) {
                  vel = cct+cst*sqrt(age[i][ipt])+clt*age[i][ipt];
                  T += dis[i][ipt]/vel;
               }
               vel = dist[i]/T;
               rsd1[icc][ics][icl][i] = fabs(vel-grv[i]);
               rsd2[icc][ics][icl][i] = pow(rsd1[icc][ics][icl][i],2);
            }
          //discard outliers
            for(ne=0,rsdt=0,i=0;i<iper;i++) if(rsd1[icc][ics][icl][i]>0) { rsdt+=rsd2[icc][ics][icl][i]; ne++; }
            rsdt /= ne;
            for(rstd=0,i=0;i<iper;i++) if(rsd1[icc][ics][icl][i]>0) rstd+=pow((rsd2[icc][ics][icl][i]-rsdt),2);
            rstd = sqrt(rstd/(ne-1));
        // cout<<"before discard: "<<ne<<endl;
            if(rsd1[icc][ics][icl][0]>0 && rsd2[icc][ics][icl][0]<rsdt+rstd*2.5) ne=1;
            else { rsd1[icc][ics][icl][0]=0.; rsd2[icc][ics][icl][0]=0.; ne=0; }
            for(i=1;i<iper;i++) if(rsd1[icc][ics][icl][i]>0 && rsd2[icc][ics][icl][i]<rsdt+rstd*2.5) { rsd1[icc][ics][icl][0]+=rsd1[icc][ics][icl][i]; rsd2[icc][ics][icl][0]+=rsd2[icc][ics][icl][i]; ne++; } //else {cout<<rsd[icc][ics][icl][i]<<" ";} cout<<endl;
            rsd1[icc][ics][icl][0] = rsd1[icc][ics][icl][0]/(ne-1);
            rsd2[icc][ics][icl][0] = sqrt(rsd2[icc][ics][icl][0]/(ne-1));
        //    rsdm = rsd[icc][ics][icl][0]/sqrt((double)ne);
      //cout<<"final: "<<rsd[icc][ics][icl][0]<<endl;
            fprintf(ff, "%f %f %f %f %f %d\n", cct, cst, clt, rsd1[icc][ics][icl][0], rsd2[icc][ics][icl][0], ne );
            if(rsdmin>rsd2[icc][ics][icl][0]) {
               rsdmin = rsd2[icc][ics][icl][0];
               *cc = cct;
               *cs = cst;
               *cl = clt;
               for(i=0;i<iper;i++) 
                  if(rsd1[icc][ics][icl][i]>0 && rsd2[icc][ics][icl][i]<rsdt+rstd*2.5) flag[i] = 1;
                  else flag[i] = 0;
//cout<<"before: "<<iper<<"  after: "<<ne<<"  rsd: "<<rsdmin<<endl;
             }
         }
      }
   }
   fclose(ff);
}

void Search_ph_age(float per, float *phv, float *lon1, float *lat1, float *lon2, float *lat2, int iper, float *cc, float *cs,  float *cl, int *flag, char *fname) {
   cout<<"Grid Searching for phase vel..."<<endl;
   float ccl, cch, csl, csh, cll, clh;
   if(per==6) { ccl=2.55, cch=2.8, csl=-0.55, csh=-0.4, cll=0., clh=0.06; }
   //if(per==6) { ccl=2.2, cch=2.2, csl=-0.297, csh=-0.297, cll=-0.0091, clh=-0.0091; }
   else if(per==6.5) { ccl=2.65, cch=2.9, csl=-0.3, csh=-0.15, cll=-0.1, clh=-0.04; }
   else if(per==7.0) { ccl=2.9, cch=3.15, csl=-0.15, csh=0., cll=-0.16, clh=-0.1; }
   else if(per==7.5) { ccl=3.1, cch=3.6, csl=-0.35, csh=-0.05, cll=-0.2, clh=0.1; }
   else if(per==8.0) { ccl=3.2, cch=3.45, csl=0.1, csh=0.25, cll=-0.18, clh=-0.12; }
   else if(per==8.5) { ccl=3.1, cch=3.6, csl=-0.05, csh=0.25, cll=-0.25, clh=0.05; }
   else if(per==9.0) { ccl=3.1, cch=3.6, csl=0., csh=0.3, cll=-0.2, clh=0.1; }
   else if(per==10.0) { ccl=3.2, cch=3.45, csl=0.3, csh=0.45, cll=-0.13, clh=-0.07; }
   else if(per==11.0) { ccl=3.1, cch=3.6, csl=0.3, csh=0.6, cll=-0.3, clh=0.; }
   else if(per==12.0) { ccl=3.13, cch=3.38, csl=0.5, csh=0.65, cll=-0.17, clh=-0.11; }
   else if(per==13.0) { ccl=3.3, cch=3.55, csl=0.25, csh=0.4, cll=-0.1, clh=-0.04; }
   else if(per==14.0) { ccl=3.3, cch=3.55, csl=0.2, csh=0.35, cll=-0.07, clh=-0.01; }
   else if(per==15.0) { ccl=3.35, cch=3.6, csl=0.2, csh=0.35, cll=-0.07, clh=-0.01; }
   else if(per==16.0) { ccl=3.45, cch=3.7, csl=0., csh=0.15, cll=0.01, clh=0.07; }
   else if(per==17.0) { ccl=3.3, cch=3.8, csl=-0.05, csh=0.25, cll=-0.1, clh=0.2; }
   else if(per==18.0) { ccl=3.4, cch=3.65, csl=0.13, csh=0.28, cll=-0.05, clh=0.01; }
   else if(per==19.0) { ccl=3.2, cch=3.7, csl=0.2, csh=0.5, cll=-0.2, clh=0.1; }
   else { ccl=3.2, cch=3.7, csl=0.25, csh=0.55, cll=-0.2, clh=0.1; }

   float cct, cst, clt;
   int i, ipt, ne, icc, ics, icl, nstep = 100, npt[iper];
   float *age[iper], *dis[iper], dist[iper];
   for(i=0;i<iper;i++) {age[i] = (float *) malloc (1000 * sizeof(float)); dis[i] = (float *) malloc (1000 * sizeof(float));}
//extract dist and age info for each path
   for(i=0;i<iper;i++) get_path_age(lon1[i], lat1[i], lon2[i], lat2[i], age[i], dis[i], &npt[i], &dist[i], fname);

   float rsd1[nstep][nstep][nstep][iper], rsd2[nstep][nstep][nstep][iper], rsdmin=100., rsdt, rsdm, rstd, T, vel;
   char name[100];
   sprintf(name,"rsdp_%.1f.txt",per);
   FILE *ff = fopen(name,"w");
   for(icc=0;icc<nstep;icc++) {
cout<<icc<<"/"<<nstep<<endl;
      cct = ccl + (cch-ccl)/(nstep-1)*icc;
      for(ics=0;ics<nstep;ics++) {
         cst = csl + (csh-csl)/(nstep-1)*ics;
         for(icl=0;icl<nstep;icl++) {
            clt = cll + (clh-cll)/(nstep-1)*icl;
          //calc rsd
            for(i=0;i<iper;i++) {
               if(npt[i]==0) { rsd1[icc][ics][icl][i] = -1.; continue; }
               T=0.;
               for(ipt=0;ipt<npt[i];ipt++) {
                  vel = cct+cst*sqrt(age[i][ipt])+clt*age[i][ipt];
                  T += dis[i][ipt]/vel;
               }
               vel = dist[i]/T;
               rsd1[icc][ics][icl][i] = fabs(vel-phv[i]);
               rsd2[icc][ics][icl][i] = pow(rsd1[icc][ics][icl][i],2);
            }
          //discard outliers
            for(ne=0,rsdt=0,i=0;i<iper;i++) if(rsd1[icc][ics][icl][i]>0) { rsdt+=rsd2[icc][ics][icl][i]; ne++; }
            rsdt /= ne;
            for(rstd=0,i=0;i<iper;i++) if(rsd1[icc][ics][icl][i]>0) rstd+=pow((rsd2[icc][ics][icl][i]-rsdt),2);
            rstd = sqrt(rstd/(ne-1));
        // cout<<"before discard: "<<ne<<endl;
            if(rsd1[icc][ics][icl][0]>0 && rsd2[icc][ics][icl][0]<rsdt+rstd*2.5) ne=1;
            else { rsd1[icc][ics][icl][0]=0.; rsd2[icc][ics][icl][0]=0.; ne=0; }
            for(i=1;i<iper;i++) if(rsd1[icc][ics][icl][i]>0 && rsd2[icc][ics][icl][i]<rsdt+rstd*2.5) { rsd1[icc][ics][icl][0]+=rsd1[icc][ics][icl][i]; rsd2[icc][ics][icl][0]+=rsd2[icc][ics][icl][i]; ne++; } //else {cout<<rsd[icc][ics][icl][i]<<" ";} cout<<endl;
            rsd1[icc][ics][icl][0] = rsd1[icc][ics][icl][0]/(ne-1);
            rsd2[icc][ics][icl][0] = sqrt(rsd2[icc][ics][icl][0]/(ne-1));
           //if(fabs(rsd1[icc][ics][icl][0]-rsd2[icc][ics][icl][0])/rsd2[icc][ics][icl][0]>0.5) {
              //for(i=1;i<iper;i++) {
              //   if(rsd1[icc][ics][icl][i]>0 && rsd2[icc][ics][icl][i]<rsdt+rstd*2.5) {
              //       cout<<"L1: "<<rsd1[icc][ics][icl][i]<<"  L2: "<<rsd2[icc][ics][icl][i]<<endl;}
              //   }
              //cout<<"L1 final: "<<rsd1[icc][ics][icl][0]<<"  L2 final: "<<rsd2[icc][ics][icl][0]<<endl; 
              
          // }
        //    rsdm = rsd[icc][ics][icl][0]/sqrt((double)ne);
      //cout<<"final: "<<rsd[icc][ics][icl][0]<<endl;
            fprintf(ff, "%f %f %f %f %f %d\n", cct, cst, clt, rsd1[icc][ics][icl][0], rsd2[icc][ics][icl][0], ne );
            if(rsdmin>rsd2[icc][ics][icl][0]) {
               rsdmin = rsd2[icc][ics][icl][0];
               *cc = cct;
               *cs = cst;
               *cl = clt;
            for(i=0;i<iper;i++)
                  if(rsd1[icc][ics][icl][i]>0 && rsd2[icc][ics][icl][i]<rsdt+rstd*2.5) flag[i] = 1;
                  else flag[i] = 0;
//cout<<"before: "<<iper<<"  after: "<<ne<<"  rsd: "<<rsdmin<<endl;
             }
         }
      }
   }
   fclose(ff);
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

int discard_list(char *sta1, char *sta2, float pero, int flag) {
   if(!strcmp(sta1,"J23A") && !strcmp(sta2,"J38A") && pero==9) return 1;
   if(!strcmp(sta1,"J23A") && !strcmp(sta2,"J54A") && pero>=13) return 1;
   if(!strcmp(sta1,"J28A") && !strcmp(sta2,"J47A") && pero>=15) return 1;
   if(!strcmp(sta1,"J28A") && !strcmp(sta2,"J68A") && pero<=10) return 1;
   if(!strcmp(sta1,"J29A") && !strcmp(sta2,"J47A") && pero==18) return 1;
   if(!strcmp(sta1,"J29A") && !strcmp(sta2,"J63A") && (pero>=18 || pero==7.5)) return 1;
   if(!strcmp(sta1,"J30A") && !strcmp(sta2,"J47A") && pero==14) return 1;
   if(!strcmp(sta1,"J31A") && !strcmp(sta2,"J54A")) return 1;
   if(!strcmp(sta1,"J38A") && !strcmp(sta2,"J63A") && pero<=7) return 1;
   if(!strcmp(sta1,"J38A") && !strcmp(sta2,"J67A") && flag<0) return 1;
   if(!strcmp(sta1,"J45A") && !strcmp(sta2,"J67A") && pero==14) return 1;
   if(!strcmp(sta1,"J45A") && !strcmp(sta2,"J55A")) return 1;
   if(!strcmp(sta1,"J46A") && !strcmp(sta2,"J63A") && pero==15) return 1;
   if(!strcmp(sta1,"J52A") && !strcmp(sta2,"J63A") && ( pero==12 || pero==13 )) return 1;
//cout<<"Not discarded: "<<sta1<<" "<<sta2<<" "<<pero<<" "<<flag<<endl;
   return 0;
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
   float agemn, agemx, ageav, stdtmp, ageavgp[NPTH], ageavgg[NPTH], stdp[NPTH], stdg[NPTH], snr[NPTH], distfact[NPTH], azii[NPTH];
   double dist, azi;
   float  pero=atof(arg[4]);

   cout<<"Period: "<<pero<<endl;
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
      //if(azi<45 || azi>135&&azi<225 || azi>315) continue; // control the path direction here
      azii[iper] = azi;
      cout<<"Extracting phv and grv of path "<<sta1[ipth]<<" - "<<sta2[ipth]<<endl;
      get_snr(snrpf[ipth], pero, &snrpos);
      get_snr(snrnf[ipth], pero, &snrneg);
      if(snrpos < 5 && snrneg < 5) continue;
      if(snrpos > 5 && snrneg > 5 && ph_match(dispnf[ipth], disppf[ipth], pero)) {
            snr[iper] = (snrpos+snrneg)/1.41;
            get_vel(dispsf[ipth], pero, &grv[iper], &phv[iper]);
            if(grv[iper]<0) continue;
            flag = 0;
      }
      else if(snrpos > snrneg) { get_vel(disppf[ipth], pero, &grv[iper], &phv[iper]); snr[iper] = snrpos; flag = 1;}
      else { get_vel(dispnf[ipth], pero, &grv[iper], &phv[iper]); snr[iper] = snrneg; flag = -1; }
      if(discard_list(sta1[ipth], sta2[ipth], pero, flag)) continue;
    //cout<<pero<<" "<<phv[iper]<<" "<<grv[iper]<<" "<<flag<<" "<<snrpos<<" "<<snrneg<<endl;
      iper++;
   }
  
cout<<"npath: "<<iper<<endl;
   float ccp, csp, clp;
   float ccg, csg, clg;
   int flagp[iper], flagg[iper];
   Search_ph_age(pero, phv, lon1, lat1, lon2, lat2, iper, &ccp, &csp, &clp, flagp, arg[3]);
   Search_gr_age(pero, grv, lon1, lat1, lon2, lat2, iper, &ccg, &csg, &clg, flagg, arg[3]);
exit(0);
   for(i=0;i<iper;i++) calc_age(lon1[i], lat1[i], lon2[i], lat2[i], arg[3], ccp, csp, clp, &ageavgp[i], &stdp[i]);
   for(i=0;i<iper;i++) calc_age(lon1[i], lat1[i], lon2[i], lat2[i], arg[3], ccg, csg, clg, &ageavgg[i], &stdg[i]);
   cout<<"From Phv: "<<ccp<<" "<<csp<<" "<<clp<<endl;
   cout<<"From Grv: "<<ccg<<" "<<csg<<" "<<clg<<endl;


   sprintf(buff, "Age_Grv_%.1fsec_GS", pero);
   ff=fopen(buff, "w");
   for(i=0;i<iper;i++) fprintf(ff, "%f %f %f %f %f %f %s %f %f %s %f %f %f %d\n", pero, grv[i], ageavgg[i], stdg[i], snr[i], distfact[i], stan[ista1[i]], slon[ista1[i]], slat[ista1[i]], stan[ista2[i]], slon[ista2[i]], slat[ista2[i]], azii[i], flagg[i]);
   fclose(ff);

   sprintf(buff, "Age_Phv_%.1fsec_GS", pero);
   ff=fopen(buff, "w");
   for(i=0;i<iper;i++) fprintf(ff, "%f %f %f %f %f %f %s %f %f %s %f %f %f %d\n", pero, phv[i], ageavgp[i], stdp[i], snr[i], distfact[i], stan[ista1[i]], slon[ista1[i]], slat[ista1[i]], stan[ista2[i]], slon[ista2[i]], slat[ista2[i]], azii[i], flagp[i]);
   fclose(ff);

  return 1;
}
