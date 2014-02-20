#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <cmath>
using namespace std;

#define min(a,b) (((a) > (b)) ? (a) : (b))

#define NNODE 10000

int main (int argc, char *argv[])
{
   if(argc!=3){
      cout<<"usage: "<<argv[0]<<" [lst of model maps (per mfile)] [out_dir]"<<endl;
      return 0;
   }

   FILE *ff;
   char flst[100][100], buff[300], ctmp[100];
   int i, j, nf, itmp;
   float perlst[100], ftmp;

//-----read in per-file list-----//
   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Cannot open file "<<argv[1]<<endl;
      return -1;
   }
   //-----Insertion Sort-----//
   for(nf=0;;nf++) {
      if(fgets(buff,300,ff)==NULL) break;
      sscanf(buff,"%f %s", &ftmp, ctmp);
      for(i=nf;i>0 && ftmp<perlst[i-1];i--) { 
         perlst[i]=perlst[i-1]; sscanf(flst[i-1],"%s",flst[i]); 
      }
      perlst[i]=ftmp;
      sscanf(ctmp,"%s",flst[i]);
      //cout<<perlst[nf]<<" "<<flst[nf]<<endl;
   }
   fclose(ff);
   cout<<nf<<" model maps read in from "<<argv[1]<<endl;

   float lonlst[nf][NNODE], latlst[nf][NNODE], vellst[nf][NNODE];
   float lonmn=360., lonmx=0., latmn=90., latmx=-90., dgrd=100., grd;
   int nn[nf];
//-----read in model maps-----//
   for(i=0;i<nf;i++) {
      if((ff=fopen(flst[i],"r"))==NULL) {
         cout<<"Cannot open file "<<flst[i]<<". Skipped."<<endl;
         continue;
      }
      for(j=0;;j++) {
         if(fgets(buff,300,ff)==NULL) break;
         sscanf(buff,"%f %f %f", &lonlst[i][j], &latlst[i][j], &vellst[i][j]);
         if(lonlst[i][j]<0.) lonlst[i][j] += 360.;
         if(lonlst[i][j]<lonmn) lonmn = lonlst[i][j];
         if(lonlst[i][j]>lonmx) lonmx = lonlst[i][j];
         if(latlst[i][j]<latmn) latmn = latlst[i][j];
         if(latlst[i][j]>latmx) latmx = latlst[i][j];
      }
      fclose(ff);
      nn[i] = j; itmp = (int)sqrt(nn[i]);
      for(j=0;j<nn[i]-1;j=j+itmp) {
         grd = min(fabs(lonlst[i][j+1]-lonlst[i][j]), fabs(latlst[i][j+1]-latlst[i][j]));
         if(dgrd>grd) dgrd=grd;
      }
   }
   cout<<"dnode: "<<dgrd<<"  region: "<<lonmn<<"-"<<lonmx<<" "<<latmn<<"-"<<latmx<<endl;
   int nlon = (int)ceil((lonmx-lonmn)/dgrd+1);
   int nlat = (int)ceil((latmx-latmn)/dgrd+1);
   int ilon, ilat, iper;
   float vel[nlon][nlat][nf];
   for(iper=0;iper<nf;iper++) {
      for(i=0;i<nlon;i++) for(j=0;j<nlat;j++) vel[i][j][iper] = -1.;
      for(i=0;i<nn[iper];i++) {
         ilon = (int)((lonlst[iper][i]-lonmn)/dgrd+0.5);
         ilat = (int)((latlst[iper][i]-latmn)/dgrd+0.5);
         vel[ilon][ilat][iper] = vellst[iper][i];
      }
   }

   float flon, flat;
   sprintf(buff,"mkdir -p %s",argv[2]);
   system(buff);
   for(ilon=0;ilon<nlon;ilon++) {
      flon = ilon*dgrd+lonmn;
      for(ilat=0;ilat<nlat;ilat++) {
         flat = ilat*dgrd+latmn;
         sprintf(ctmp,"%s/%.2f_%.2f.disp",argv[2],flon,flat);
         for(i=0,iper=0;iper<nf;iper++) if(vel[ilon][ilat][iper]>0){i=1; break;}
         if(i==0) continue;
         ff=fopen(ctmp,"w");
         for(iper=0;iper<nf;iper++) {
            if(vel[ilon][ilat][iper]<0.) continue;
            fprintf(ff, "%f %f\n", perlst[iper], vel[ilon][ilat][iper]);
         }
         fclose(ff);
      }
   }
 
   return 1;
}

