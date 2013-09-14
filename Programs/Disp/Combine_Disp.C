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

main(int na, char *arg[])
{
  if (na!=4)
    {
      cout<<"usage: "<<arg[0]<<" [station.lst] [sta1-sta2-Dispf-SNRf list] [out_name]"<<endl;
      return 0;
    }
   
   FILE *ff, *fsnr, *fdsp;
   char buff[300], ctmp[100];
   char stan[NSTA][6], sta1[NPTH][6], sta2[NPTH][6], dispf[NPTH][100], snrf[NPTH][100];
   int i, j, iper, ipth, npth, nsta, itmp;
   float slon[NSTA], slat[NSTA], ftmp, snrpos, snrneg;
   float lon1, lon2, lat1, lat2;
   float grv[NPTH*20], phv[NPTH*20], per[NPTH*20], pertmp;
   double dist;

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
      lon1 = slon[i]; lat1 = slat[i];
      for(i=0;i<nsta;i++) if(strcmp(sta2[ipth],stan[i])==0) break;
      if(i==nsta) continue;
      lon2 = slon[i]; lat2 = slat[i];
      calc_dist(lat1, lon1, lat2, lon2, &dist);
      if((fsnr=fopen(snrf[ipth],"r"))==NULL) {
         cout<<"Cannot open file "<<snrf[ipth]<<endl;
         continue;
      }
      if((fdsp=fopen(dispf[ipth],"r"))==NULL) {
         cout<<"Cannot open file "<<dispf[ipth]<<endl;
         continue;
      }
      for(i=0;;i++) {
         if(fgets(buff, 300, fsnr)==NULL) break;
         if(sscanf(buff,"%f %f %f %f %f", &pertmp, &ftmp, &snrpos, &ftmp, &snrneg)!=5) {cout<<"Wrong format in file "<<snrf[ipth]<<endl; return -1;}
         if(fgets(buff, 300, fdsp)==NULL) break;
         if(sscanf(buff,"%d %f %f %f %f", &itmp, &ftmp, &per[iper], &grv[iper], &phv[iper])!=5) {cout<<"Wrong format in file "<<dispf[ipth]<<endl; return -1;}
         if(pertmp!=per[iper]) {cout<<"Mismatch between disp and snr file"<<endl; break;}
         if(per[iper]*10. > dist) continue;
         if( snrpos+snrneg < 12 ) continue;
   if(per[iper]<15 && grv[iper]>4) cout<<snrf[ipth]<<endl;
         iper++;
      }
      fclose(fdsp); fclose(fsnr);
   }
cout<<"nper: "<<iper<<endl;
  
   sprintf(buff, "%s.gr.disp", arg[3]);
   ff=fopen(buff, "w");
   for(i=0;i<iper;i++) fprintf(ff, "%f %f\n", per[i], grv[i]);
   fclose(ff);
   sprintf(buff, "%s.ph.disp", arg[3]);
   ff=fopen(buff, "w");
   for(i=0;i<iper;i++) fprintf(ff, "%f %f\n", per[i], phv[i]);
   fclose(ff);

  return 1;
}
