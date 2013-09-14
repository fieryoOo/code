#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
using namespace std;

int main(int argc, char *argv[])
{
   if(argc!=3){
      cout<<"Usage: Map_Correlation [in_HD_file_1] [in_HD_file_2]"<<endl;
      return 0;
   }

   FILE *ff;
   char buff[300];
   float longmin=360, longmax=0, latimin=90, latimax=-90, tmp;
   float dgr[2], grdlon[200000][2], grdlat[200000][2];
   double grddat[200000][2], corr=0, weight1=0, weight2=0;
   int i, j, ii, jj, ngrd[2];

   for(i=0;i<2;i++){
      if((ff=fopen(argv[i+1],"r")) == NULL){
         cout<<"Can't open file "<<argv[i+1]<<" for read!"<<endl;
         return 0;
      }
      for(ngrd[i]=0;;ngrd[i]++){
         if(fgets(buff, 300, ff) == NULL) break;
         sscanf(buff,"%f %f %lf", &grdlon[ngrd[i]][i], &grdlat[ngrd[i]][i], &grddat[ngrd[i]][i]);
         if(longmin>grdlon[ngrd[i]][i]) longmin=grdlon[ngrd[i]][i];
         if(longmax<grdlon[ngrd[i]][i]) longmax=grdlon[ngrd[i]][i];
         if(latimin>grdlat[ngrd[i]][i]) latimin=grdlat[ngrd[i]][i];
         if(latimax<grdlat[ngrd[i]][i]) latimax=grdlat[ngrd[i]][i];
      }
      fclose(ff);
      if(ngrd[i]<100){
         cout<<"Only "<<ngrd[i]<<" points in file "<<argv[i+1]<<endl;
         return 0;
      }

      dgr[i]=100;
      for(j=0;j<99;j++){
         tmp=fabs(grdlon[j][i]-grdlon[j+1][i]);
         if(tmp==0) tmp=fabs(grdlat[j][i]-grdlat[j+1][i]);
         if(tmp==0) continue;
         if(dgr[i]>tmp) dgr[i]=tmp;
      }
      if(dgr[i]==100){
         cout<<"Wrong HD file!"<<endl;
         return 0;
      }
   }

   if(dgr[0]!=dgr[1])
      cout<<"Grid res do not match. The smaller is used automatically"<<endl;
   if(dgr[0]>dgr[1]) dgr[0]=dgr[1];

//   longmin=floor(longmin/dgr[0])*dgr[0];
//   latimin=floor(latimin/dgr[0])*dgr[0];
//   longmax=ceil(longmax/dgr[0])*dgr[0];
//   latimax=ceil(latimax/dgr[0])*dgr[0];
   int npts_long=int((longmax-longmin)/dgr[0]+1);
   int npts_lati=int((latimax-latimin)/dgr[0]+1);
   int flag[npts_long][npts_lati];
   double dat[npts_long][npts_lati][2];
   double dat_avg[2];

   for(i=0;i<npts_long;i++) for(j=0;j<npts_lati;j++) flag[i][j] = -1;
   for(i=0;i<2;i++) {
      dat_avg[i] = 0.;
      for(j=0;j<ngrd[i];j++){
         dat_avg[i] += grddat[j][i];
         ii=int((grdlon[j][i]-longmin)/dgr[0]+0.5);
         jj=int((grdlat[j][i]-latimin)/dgr[0]+0.5);
         dat[ii][jj][i]=grddat[j][i];
         flag[ii][jj]+=1;
      }
      dat_avg[i] /= ngrd[i];
   }

   for(i=0;i<npts_long;i++) for(j=0;j<npts_lati;j++) {
      if(flag[i][j]<1) continue;
      if(flag[i][j]>1){
         cout<<"Incorrect grid size!"<<endl;
         return 0;
      }
      corr+=(dat[i][j][0]-dat_avg[0])*(dat[i][j][1]-dat_avg[1]);
      weight1+=pow(dat[i][j][0]-dat_avg[0],2);
      weight2+=pow(dat[i][j][1]-dat_avg[1],2);
   }
   corr=corr/sqrt(weight1*weight2);
   cout<<corr<<endl;
   return 1;
}
