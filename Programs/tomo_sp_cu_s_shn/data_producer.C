#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <cmath>
using namespace std;

#define NSTA 3000
#define NOSTA 1000

int main (int argc, char *argv[])
{
   if(argc!=3){
      cout<<"usage: "<<argv[0]<<" [lst of csta_gr_ph_amp_map files] [outf_label]"<<endl;
      return 0;
   }
   FILE *flst, *ff;
   int i, j, k, ic, nc, nfile, npath, ncsta, nsta[NSTA];
   float clon[NSTA], clat[NSTA], olon[NSTA][NOSTA], olat[NSTA][NOSTA], gvel[NSTA][NOSTA], pvel[NSTA][NOSTA], dtmp;
   char csta[NSTA][6], osta[NSTA][NOSTA][6], file[100], buff[300];
  //-----read in all csta-osta pairs and discard(warn) duplicate paths-----//
   if((flst = fopen(argv[1],"r")) == NULL) {
      cout<<"Cannot open file "<<argv[1]<<endl;
      return 0;
   }
   for(nc=0,nfile=0,npath=0;;) {
      if((fgets(buff, 300, flst)) == NULL) break;
      sscanf(buff,"%s",file);
      if((ff = fopen(file,"r")) == NULL) {
         cout<<"Cannot open file "<<file<<". Skipped."<<endl;
         continue;
      }
      if((fgets(buff, 300, ff))==NULL) continue;
      if((sscanf(buff, "%f %f %f %f %f %f %s", &clon[nc], &clat[nc], &dtmp, &dtmp, &dtmp, &dtmp, &csta[nc]))!=7) {
         cout<<"Wrong format in file "<<file<<". Skipped."<<endl;
         continue;;
      }
      nfile++;
      for(i=nc-1;i>=0;i--) if(strcmp(csta[i], csta[nc]) == 0) break;
      if(i>=0) ic=i;
      else { ic=nc; nsta[ic]=0; nc++; }
      for(i=nsta[ic];;) {
         if((fgets(buff, 300, ff))==NULL) break;
         if((sscanf(buff, "%f %f %f %f %f %f %s", &olon[ic][i], &olat[ic][i], &dtmp, &gvel[ic][i], &pvel[ic][i], &dtmp, &osta[ic][i]))!=7) continue;
         for(j=i-1;j>=0;j--) if(strcmp(osta[ic][i],osta[ic][j])==0) break;
         if(j>=0) {
            if(fabs(pvel[ic][i]-pvel[ic][j])>0.01)cout<<"Warning: over 1 different phvels found for path "<<csta[ic]<<" - "<<osta[ic][i]<<endl;
            continue;
         }
         for(j=nc-1;j>=0;j--) {
            if(j==ic) continue;
            if(strcmp(osta[ic][i],csta[j])==0) break;
         }
         if(j>=0) {
            for(k=nsta[j]-1;k>=0;k--) if(strcmp(csta[ic],osta[j][k])==0) break;
            if(k>=0) {
               if(fabs(pvel[ic][i]-pvel[j][k])>0.01)cout<<"Warning: over 1 different phvels found for path "<<csta[ic]<<" - "<<osta[ic][i]<<endl;
               continue;
            }
         }
         i++; npath++;
      }
      fclose(ff);
      nsta[ic] = i;
   }
   fclose(flst);
   cout<<nfile<<" files read in"<<endl;
   cout<<"# of center stations: "<<nc<<endl;
   cout<<"# of paths: "<<npath<<endl;
  //-----output data file-----//
   sprintf(buff, "%s_gr", argv[2]);
   flst = fopen(buff,"w");
   for(ic=0,j=1;ic<nc;ic++) for(i=0;i<nsta[ic];i++,j++)
      fprintf(flst, "%6d %10.5f %10.5f %10.5f %10.5f %8.5f   1 %6s %6s\n", j, clat[ic], clon[ic], olat[ic][i], olon[ic][i], gvel[ic][i], csta[ic], osta[ic][i]);
   fclose(flst);
   sprintf(buff, "%s_ph", argv[2]);
   flst = fopen(buff,"w");
   for(ic=0,j=1;ic<nc;ic++) for(i=0;i<nsta[ic];i++,j++)
      fprintf(flst, "%6d %10.5f %10.5f %10.5f %10.5f %8.5f   1 %6s %6s\n", j, clat[ic], clon[ic], olat[ic][i], olon[ic][i], pvel[ic][i], csta[ic], osta[ic][i]);
   fclose(flst);

   return 1;
}

