#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <time.h>
using namespace std;

int calc_dist(double lati1, double long1, double lati2, double long2, double *dist)
{
  int i,flag=0;
  double latio1,longo1,dlt_lati,dlt_long;
  double Ra,Rb,R,Rx,f;
  double U1,U2;
  double Ds,alpha2;
  double ctr_angl;
  double pi;
  double cv,cv1,cv2,cv3,cv4,cv5,cvC,numda1;
  double numda;
  double mius, cvA, cvB, deltacv;

  pi=4.0*atan(1.0);
  Ra = 6378.137; Rb = 6356.7523142;
  f = 1/298.257223563;
begin:
  long1=long1-(int)floor(long1)/360*360;
  if(long1<0) long1+=360;
  long2=long2-(int)floor(long2)/360*360;
  if(long2<0) long2+=360;
  if(lati1==-lati2 && fabs(long1-long2)==180){
  *dist=20003.917357;
  return 1;
  }
  dlt_lati=fabs(lati2-lati1);
  dlt_long=long2-long1;

if (dlt_long > 180.000)  dlt_long = 360.000000 - dlt_long;
if (dlt_long < -180.000) dlt_long = 360.000 - fabs(dlt_long);
dlt_long = fabs(dlt_long);

U1 = atan((1-f)*tan(lati1/180*pi));
U2 = atan((1-f)*tan(lati2/180*pi));
dlt_long = dlt_long*pi/180;

numda = dlt_long;
numda1 = numda;
i=0;
do {
i++;
numda = numda1;
cv1 =  sqrt( (cos(U2)*sin(numda))*(cos(U2)*sin(numda))+ (cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda))*(cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda)) );
cv2 = sin(U1)*sin(U2)+ cos(U1)*cos(U2)*cos(numda);
cv = atan2(cv1,cv2);
if(cv==0)cv+=0.00000000001;
cv3 = cos(U1)*cos(U2)*sin(numda)/sin(cv);
cv4 = 1 - cv3*cv3;
if(cv4==0)cv4+=0.00000000001;
cv5 = cos(cv) - 2*sin(U1)*sin(U2)/cv4;
cvC = f/16*cv4*(4 + f*(4 - 3*cv4));
numda1 = dlt_long + (1-cvC)*f*cv3*(cv + cvC*cv1*(cv5 + cvC*cv2*(-1 +2*cv5*cv5)));
if(i>50){
flag=1;
latio1=lati1*pi/180.;
Rx=sqrt((pow(Ra*Ra*cos(latio1),2)+pow(Rb*Rb*sin(latio1),2))/(pow(Ra*cos(latio1),2)+pow(Rb*sin(latio1),2)));
Rx=pi*(Rx+Ra)/2.;
latio1=lati1;
longo1=long1;
lati1=lati2;
long1=long2;
lati2=-latio1;
long2=longo1+180;
goto begin;
}
} while (fabs(numda - numda1) > 1e-10);
mius = cv4*(Ra*Ra - Rb*Rb)/(Rb*Rb);
cvA = 1+mius/16384*(4096 + mius*(-768 + mius*(320 - 175*mius)));
cvB = mius/1024*(256+ mius*(-128 + mius*(74 - 47*mius)));
deltacv = cvB*cv1*(cv5 +cvB/4*(cv2*(-1 + 2*cv5*cv5)-cvB/6*cv5*(-3+4*cv1*cv1)*(-3+4*cv5*cv5) ));
*dist = Rb * cvA *(cv - deltacv);

if(flag==1){
alpha2=atan2(cos(U1)*sin(dlt_long), -sin(U1)*cos(U2) + cos(U1)*sin(U2)*cos(dlt_long))*180/pi;
if( fabs(long2-long1)>180 ) { alpha2 = 360-alpha2;}
if( long2 < long1 ) { alpha2 = 360-alpha2; }
if(alpha2>180) alpha2=360-alpha2;
alpha2=fabs(90-alpha2);
Ds=Rx*(90-alpha2)/90.+(Ra+Rb)*pi*alpha2/180.;
*dist=Ds-*dist;
}

return 1;
}


int calc_azimuth(double lati1, double long1, double lati2, double long2, double *alpha1)
{
  double dlt_lati,dlt_long;
  double Ra,Rb,f;
  double U1,U2;
  double pi;

  pi=4.0*atan(1.0);
  Ra = 6378.137; Rb = 6356.7523142;
  f = 1/298.257223563;
  long1=long1-(int)floor(long1)/360*360;
  if(long1<0) long1+=360;
  long2=long2-(int)floor(long2)/360*360;
  if(long2<0) long2+=360;
  if(lati1==-lati2 && fabs(long1-long2)==180){
  *alpha1=999;
  return 1;
  }
  dlt_lati=fabs(lati2-lati1);
  dlt_long=long2-long1;

if (dlt_long > 180.000)  dlt_long = 360.000000 - dlt_long;
if (dlt_long < -180.000) dlt_long = 360.000 - fabs(dlt_long);
dlt_long = fabs(dlt_long);

U1 = atan((1-f)*tan(lati1/180*pi));
U2 = atan((1-f)*tan(lati2/180*pi));
dlt_long = dlt_long*pi/180;

*alpha1=atan2(cos(U2)*sin(dlt_long), cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(dlt_long))*180/pi;
//*alpha2=atan2(cos(U1)*sin(dlt_long), -sin(U1)*cos(U2) + cos(U1)*sin(U2)*cos(dlt_long))*180/pi;
if( fabs(long2-long1)>180 ) { *alpha1 = 360-*alpha1; }
if( long2 < long1 ) { *alpha1 = 360-*alpha1; }

return 1;
}

int Check_azi_cov(double *lon, double *lat, int nsta, double mdis, double *flag)
{
   int i, j, jj, k, ndata, ntmp;
   double azi[300], azimin, azitmp, dist;
   double mazi=110.;

   for(i=0;i<nsta;i++){
      if(flag[i]==-1) continue;
      ndata=0;
      for(j=0;j<nsta;j++){
         if(j==i || flag[j]==-1) continue;
         calc_dist(lat[i], lon[i], lat[j], lon[j], &dist);
         if(dist>mdis) continue;
         calc_azimuth(lat[i], lon[i], lat[j], lon[j], &azi[ndata]);
         ndata++;
      }
      //cout<<"Check azi cov: "<<ndata<<" nearby stations"<<endl;
      for(j=0;j<ndata;j++) {
         azimin=azi[j]; jj=j;
         for(k=j+1;k<ndata;k++)
            if(azimin>azi[k]) {
              azimin=azi[k];
              jj=k;
            }
         if(jj==j) continue;
         azitmp=azi[j];
         azi[j]=azi[jj];
         azi[jj]=azitmp;
      }
      if(azi[0]+360.-azi[ndata-1]>mazi) {
         flag[i]=-2;
         continue;
      }
      for(j=1;j<ndata;j++)
         if(azi[j]-azi[j-1]>mazi) {
            flag[i]=-2;
            break;
         }
   }
   for(i=0;i<nsta;i++) if(flag[i]==-2) flag[i]=-1;
   return 1;
}

int calc_lplc(char *file, double *slong, double *slati, double *grdt, double *azi, double *lplc, int *nsta, double tension, double bdis)
{
   FILE *ff;
   char buff[300], buff2[300];
   float longmin=360, longmax=0, latimin=90, latimax=-90, longtmp, latitmp;
   double radius=6371.1391285, pi=4.0*atan(1.0), dgr=0.05, dattmp;
   double grdtx, grdty, dis, alpha;
   double sdat[2000], weight;
   int i, j, ii, jj, ista, nrand, ntmpp, nstl;

   if((ff=fopen(file,"r"))==NULL){
     printf("Can't open %s to read!\n",file);
     return 0;
   }
   for(*nsta=0;;){
      if(fgets(buff, 300, ff) == NULL ) break;
      if((sscanf(buff,"%lf %lf %lf", &slong[*nsta], &slati[*nsta], &sdat[*nsta]))!=3) {
         cout<<"Wrong format! Stopped: "<<buff<<endl;
         break;
      }
      if(longmin>slong[*nsta]) longmin=slong[*nsta];
      if(longmax<slong[*nsta]) longmax=slong[*nsta];
      if(latimin>slati[*nsta]) latimin=slati[*nsta];
      if(latimax<slati[*nsta]) latimax=slati[*nsta];
      *nsta=*nsta+1;
   }
   fclose(ff);

   for(i=0;i<*nsta;i++) azi[i] = 0;
   Check_azi_cov(&slong[0], &slati[0], *nsta, 200, &azi[0]);
   //nstl=0;
   //for(i=0;i<*nsta;i++) if(azi[i]!=-1) nstl++;
   //cout<<"lplc: found "<<nstl<<" out of "<<*nsta<<" stations with good coverage."<<endl;

   longmin=floor(longmin/dgr)*dgr-dgr;
   latimin=floor(latimin/dgr)*dgr-dgr;
   longmax=ceil(longmax/dgr)*dgr+dgr;
   latimax=ceil(latimax/dgr)*dgr+dgr;
   int npts_long=int((longmax-longmin)/dgr+1);
   int npts_lati=int((latimax-latimin)/dgr+1);
   double dislong[npts_lati], dislati, distmp;
   double gradtr[npts_long][npts_lati], grdtmp, gradx[npts_long][npts_lati], grady[npts_long][npts_lati], grdx[*nsta], grdy[*nsta], dat[npts_long][npts_lati];

   for(i=0;i<npts_long;i++) for(j=0;j<npts_lati;j++){
      gradtr[i][j] = 1e10;
      dat[i][j] = -1;
   }
/*
   for(i=1;i<npts_lati-1;i++){
      dislati=atan(0.993277*tan((latimin+i*dgr)/180*pi))*180/pi;
      dislong[i]=radius*sin((90-dislati)/180*pi)*dgr/180*pi;
   }
   dislati=radius*dgr/180*pi;
*/
   dislati=dgr; 
   for(i=1;i<npts_lati-1;i++) dislong[i]=dgr;
   nrand=rand();
   sprintf(buff,"region_%d",nrand);
   ff=fopen(buff,"w");
   fprintf(ff, "-R%d/%d/%d/%d\0", (int)floor(longmin), (int)ceil(longmax), (int)floor(latimin), (int)ceil(latimax));
   fclose(ff);
   sprintf(buff,"/media/WORK/tianye/Test_tension/test_sinAoA/C_plot_travel %s region_%d %f %.2f %.0f", file, nrand, dgr, 0., bdis);
   system(buff);

   sprintf(buff,"%s_ts%.2f_bs%.0f.HD",file,0.,bdis);
   if((ff=fopen(buff,"r"))==NULL){
      printf("Can't open %s to read!\n",buff);
      return 0;
   }
   for(;;){
      if( fgets(buff, 300, ff) == NULL ) break;
      sscanf(buff,"%f %f %lf", &longtmp, &latitmp, &dattmp);
      i=int((longtmp-longmin)/dgr+0.1);
      j=int((latitmp-latimin)/dgr+0.1);
      if(i<0||i>=npts_long||j<0||j>=npts_lati) continue;
      dat[i][j]=dattmp;
   }
   fclose(ff);
   for(ista=0;ista<*nsta;ista++){
      if(azi[ista]==-1) continue;
      ii=int((slong[ista]-longmin)/dgr);
      jj=int((slati[ista]-latimin)/dgr);
      for(i=ii;i<ii+2;i++)
         for(j=jj;j<jj+2;j++){
            if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
            if(gradtr[i][j] != 1e10) continue;
            gradx[i][j]=(dat[i+1][j]-dat[i-1][j])/2.0/dislong[j];
            grady[i][j]=(dat[i][j+1]-dat[i][j-1])/2.0/dislati;
            gradtr[i][j]=sqrt(pow(gradx[i][j],2)+pow(grady[i][j],2));
            //azig[i][j]=pi/2-atan2(grady[i][j],gradx[i][j]);
            //if(azig[i][j]<0-1e-8)azig[i][j]+=2*pi;
            distmp=sqrt(dislong[j]*dislong[j]+dislati*dislati);
            grdtx=(dat[i+1][j-1]-dat[i-1][j+1])/2.0/distmp;
            grdty=(dat[i+1][j+1]-dat[i-1][j-1])/2.0/distmp;
            alpha=2*atan(dislati/dislong[j]);
            alpha=atan((grdtx/grdty-cos(alpha))/sin(alpha));
            grdtmp=fabs(grdty/cos(alpha));
            if(fabs((gradtr[i][j]-grdtmp)/gradtr[i][j])>0.15) gradtr[i][j]=1e20;
//            else gradtr[i][j]=(gradtr[i][j]+grdtmp)/2.;
         }
   }

   nstl=0;
   for(ista=0;ista<*nsta;ista++){
      if(azi[ista]==-1) continue;
      ii=int((slong[ista]-longmin)/dgr);
      jj=int((slati[ista]-latimin)/dgr);
      ntmpp=0; weight=0; grdx[ista]=0, grdy[ista]=0, grdt[ista]=0;
      for(i=ii;i<ii+2;i++)
         for(j=jj;j<jj+2;j++){
            if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
            if(gradtr[i][j]>=1e10) continue;
            ntmpp++;
            dis=sqrt(pow((longmin+i*dgr-slong[ista])*dislong[j],2)+pow((latimin+j*dgr-slati[ista])*dislati,2))+1e-20;
            grdx[ista]+=gradx[i][j]/dis;
            grdy[ista]+=grady[i][j]/dis;
            //grdt[ista]+=gradtr[i][j]/dis;
            //azi[ista]+=azig[i][j]/dis;
            weight+=1./dis;
         }
      if(ntmpp<2){ azi[ista]=-1; continue; }
      grdx[ista]=grdx[ista]/weight;
      grdy[ista]=grdy[ista]/weight;
      grdt[ista]=sqrt(pow(grdx[ista],2)+pow(grdy[ista],2));
      azi[ista]=pi/2-atan2(grdy[ista],grdx[ista]);
      if(azi[ista]<0-1e-8)azi[ista]+=2*pi;
      nstl++;
   }

   if(nstl<10) {
      cout<<nstl<<" stations left. Not enough for computing lplc map."<<endl;
      for(ista=0;ista<*nsta;ista++) azi[ista]=-1;
      sprintf(buff,"rm -f region_%d", nrand);
      system(buff);
      return 0;
   }

   sprintf(buff,"%s_grdx",file);
   if((ff=fopen(buff,"w"))==NULL){
      printf("Can't open %s to write!\n",buff);
      return 0;
   }
   for(ista=0;ista<*nsta;ista++) {
      if(azi[ista]==-1) continue;
      fprintf(ff,"%8.4f  %8.4f  %8g\n",slong[ista], slati[ista], grdx[ista]);
   }
   fclose(ff);

   sprintf(buff,"%s_grdy",file);
   if((ff=fopen(buff,"w"))==NULL){
      printf("Can't open %s to write!\n",buff);
      return 0;
   }
   for(ista=0;ista<*nsta;ista++) {
      if(azi[ista]==-1) continue;
      fprintf(ff,"%8.4f  %8.4f  %8g\n",slong[ista], slati[ista], grdy[ista]);
   }
   fclose(ff);

   sprintf(buff,"/media/WORK/tianye/Test_tension/test_sinAoA/C_plot_travel %s_grdx region_%d %f %.2f %.0f", file, nrand, dgr, tension, bdis);
   system(buff);
   sprintf(buff,"/media/WORK/tianye/Test_tension/test_sinAoA/C_plot_travel %s_grdy region_%d %f %.2f %.0f", file, nrand, dgr, tension, bdis);
   system(buff);

   for(i=0;i<npts_long;i++) for(j=0;j<npts_lati;j++){
      gradtr[i][j] = 1e10;
      gradx[i][j] = -1;
   }

   sprintf(buff,"%s_grdx_ts%.2f_bs%.0f.HD",file,tension,bdis);
   if((ff=fopen(buff,"r"))==NULL){
      printf("Can't open %s to read!\n",buff);
      return 0;
   }
   for(;;){
      if( fgets(buff, 300, ff) == NULL ) break;
      sscanf(buff,"%f %f %lf", &longtmp, &latitmp, &dattmp);
      i=int((longtmp-longmin)/dgr+0.1);
      j=int((latitmp-latimin)/dgr+0.1);
      if(i<0||i>=npts_long||j<0||j>=npts_lati) continue;
      gradx[i][j]=dattmp;
   }
   fclose(ff);

   sprintf(buff,"%s_grdy_ts%.2f_bs%.0f.HD",file,tension,bdis);
   if((ff=fopen(buff,"r"))==NULL){
      printf("Can't open %s to read!\n",buff);
      return 0;
   }
   for(;;){
      if( fgets(buff, 300, ff) == NULL ) break;
      sscanf(buff,"%f %f %lf", &longtmp, &latitmp, &dattmp);
      i=int((longtmp-longmin)/dgr+0.1);
      j=int((latitmp-latimin)/dgr+0.1);
      if(i<0||i>=npts_long||j<0||j>=npts_lati) continue;
      grady[i][j]=dattmp;
   }
   fclose(ff);

   Check_azi_cov(&slong[0], &slati[0], *nsta, 200, &azi[0]);

   for(ista=0;ista<*nsta;ista++){
      if(azi[ista]==-1) continue;
      ii=int((slong[ista]-longmin)/dgr);
      jj=int((slati[ista]-latimin)/dgr);
      for(i=ii;i<ii+2;i++)
         for(j=jj;j<jj+2;j++){
            if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
            if(gradtr[i][j] != 1e10) continue;
            grdtx=(gradx[i+1][j]-gradx[i-1][j])/2.0/dislong[j];
            grdty=(grady[i][j+1]-grady[i][j-1])/2.0/dislati;
            grdtmp=grdtx+grdty;
            alpha=atan(dislong[j]/dislati);
            gradtr[i][j]=-gradx[i-1][j+1]*sin(alpha)+grady[i-1][j+1]*cos(alpha)+grady[i][j+1]+gradx[i+1][j+1]*sin(alpha)+grady[i+1][j+1]*cos(alpha)-gradx[i-1][j]+gradx[i+1][j]-gradx[i-1][j-1]*sin(alpha)-grady[i-1][j-1]*cos(alpha)-grady[i][j-1]+gradx[i+1][j-1]*sin(alpha)-grady[i+1][j-1]*cos(alpha);
            gradtr[i][j]=gradtr[i][j]/8*(dislati+dislong[j])/dislati/dislong[j];
            if(fabs((grdtmp-gradtr[i][j])/gradtr[i][j])>0.5 && fabs((grdtmp-gradtr[i][j]))>0.0001) gradtr[i][j]=1e20;
            //else gradtr[i][j]=grdtmp;
         }
   }

   for(ista=0;ista<*nsta;ista++){
      if(azi[ista]==-1) continue;
      ii=int((slong[ista]-longmin)/dgr);
      jj=int((slati[ista]-latimin)/dgr);
      ntmpp=0; weight=0; lplc[ista]=0;
      for(i=ii;i<ii+2;i++)
         for(j=jj;j<jj+2;j++){
            if(i==0||j==0||i>npts_long-2||j>npts_lati-2) continue;
            if(gradtr[i][j]>=1e10) continue;
            ntmpp++;
            dis=sqrt(pow((longmin+i*dgr-slong[ista])*dislong[j],2)+pow((latimin+j*dgr-slati[ista])*dislati,2))+1e-20;
            lplc[ista]+=gradtr[i][j]/dis;
            weight+=1/dis;
         }
      if(ntmpp<2) { azi[ista]=-1; continue; }
      lplc[ista]=lplc[ista]/weight;
   }

   sprintf(buff,"rm -f %s_grdx %s_grdx*.HD %s_grdy %s_grdy*.HD region_%d", file, file, file, file, nrand);
   system(buff);

   return 1;
}

int main (int argc, char *argv[])
{
   if(argc!=4){
      cout<<"usage:csta_map_gradient [input_map] [tension] [block_size]"<<endl;
      return 0;
   }
   FILE *ff;
   char fname[100], *tmpc;
   double grdt[2000], thetat[2000], lplc[2000];
   double tlong[2000], tlati[2000], dist[2000];
   double dtmp, dis;
   double R=6367.5, pi=4.0*atan(1.0), tension=atof(argv[2]), bdis=atof(argv[3]);
   //double clon=atof(argv[4]), clat=atof(argv[5]);
   int i, nstt;

   if(bdis<0.1) bdis=0.1;

   calc_lplc(argv[1], tlong, tlati, grdt, thetat, lplc, &nstt, tension, bdis);
//for(i=1;i<nstt;i++) cout<<tlong[i]<<" "<<tlati[i]<<" "<<grdt[i]<<" "<<azit[i]<<" "<<lplc[i]<<endl;
//   if(nsta != nstt){
//      cout<<"Different station numbers from phase and amplitude map: ";
//      cout<<nsta<<" "<<nstt<<endl;
//      return 0;
//   }

   for(i=0;i<nstt;i++){
      if(thetat[i]==-1) continue;
      dtmp=pi/2-atan2(tlati[i],tlong[i]);
      if(dtmp<0-1e-8) dtmp += 2*pi;
      if(fabs(dtmp-thetat[i])>pi/2.) {
         grdt[i] = -grdt[i];
         thetat[i] -= pi;
         if(thetat[i]<0) thetat[i] += 2*pi;
      }
   }

   sprintf(fname,"%s_ts%.2f_bs%.0f_grdt",argv[1],0.,bdis);
   ff=fopen(fname,"w");
   for(i=0;i<nstt;i++){
      if(thetat[i]==-1) continue;
      fprintf(ff,"%8.4f  %8.4f  %8g\n", tlong[i], tlati[i], grdt[i]);
   }
   fclose(ff);

   sprintf(fname,"%s_ts%.2f_bs%.0f_thetat",argv[1],0.,bdis);
   ff=fopen(fname,"w");
   for(i=0;i<nstt;i++){
      if(thetat[i]==-1) continue;
      fprintf(ff,"%8.4f  %8.4f  %8g  %8g\n", tlong[i], tlati[i], thetat[i]*180./3.1415927, grdt[i]);
   }
   fclose(ff);

   sprintf(fname,"%s_ts%.2f_bs%.0f_lplc",argv[1],tension,bdis);
   ff=fopen(fname,"w");
   for(i=0;i<nstt;i++){
      if(thetat[i]==-1) continue;
      //calc_dist(tlati[i], tlong[i], clat, clon, &dis);
      //dtmp=1.4*(-0.5/grdt[i]*lplc[i]+0.5/R/tan(dis/R));
      fprintf(ff,"%8.4f  %8.4f  %8g\n", tlong[i], tlati[i], lplc[i]);
   }
   fclose(ff);

   return 1;
}
