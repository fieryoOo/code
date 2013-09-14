#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

double get_dist(double lat1,double lon1,double lat2,double lon2)
{
  double theta,pi,temp;
  double radius=6371;
  pi=4.0*atan(1.0);

  lat1=atan(0.993277*tan(lat1/180*pi))*180/pi;
  lat2=atan(0.993277*tan(lat2/180*pi))*180/pi;

  temp=sin((90-lat1)/180*pi)*cos(lon1/180*pi)*sin((90-lat2)/180*pi)*cos(lon2/180*pi)+sin((90-lat1)/180*pi)*sin(lon1/180*pi)*sin((90-lat2)/180*pi)*sin(lon2/180*pi)+cos((90-lat1)/180*pi)*cos((90-lat2)/180*pi);
  if(temp>1)
    {
      cout<<"warning cos(theta)>1 and correct to 1!!"<<temp<<endl;
      temp=1;
    }
  if(temp<-1)
    {
      cout<<"warning cos(theta)<-1 and correct to -1!!"<<temp<<endl;
      temp=-1;
    }
  theta=fabs(acos(temp));
  return theta*radius;
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

int fit_line ( double *datax, double *datay, int idata, double *slpout, double *itcptout, double *rsdout )
{
  FILE *f;
  int i,j,ii,nstep,step;
  double slp,slpmx,slpmn,slpstp,itcpt,itcptmx,itcptmn,itcptstp,slpfit,itcptfit,slpp[100],itcptt[100],rsd[100],dataxx,datayy,datamin,rsdmn,rsdfit;


  for(i=0;i<idata;i++) {
     datamin=datax[i]; ii=i;
     for(j=i+1;j<idata;j++)
        if(datamin>datax[j]){
           datamin=datax[j];
           ii=j;
        }
     if(ii==i) continue;
     dataxx=datax[i];
     datayy=datay[i];
     datax[i]=datax[ii];
     datay[i]=datay[ii];
     datax[ii]=dataxx;
     datay[ii]=datayy;
  }

  slpmx=-999999999.;slpmn=999999999.;
  itcptmx=-999999999.;itcptmn=999999999.;
  dataxx=(datax[0]+datax[1]+datax[2])/3;
  datayy=(datay[0]+datay[1]+datay[2])/3;
  for(i=idata/2;i<idata;i++){
     slp=(datay[i]-datayy)/(datax[i]-dataxx);
     //if(slp>0) continue;
     if(slpmx<slp)slpmx=slp;
     if(slpmn>slp)slpmn=slp;
    }
  dataxx=(datax[idata-1]+datax[idata-2]+datax[idata-3])/3;
  datayy=(datay[idata-1]+datay[idata-2]+datay[idata-3])/3;
  for(i=0;i<idata/2;i++){
     slp=(datayy-datay[i])/(dataxx-datax[i]);
     //if(slp>0) continue;
     if(slpmx<slp)slpmx=slp;
     if(slpmn>slp)slpmn=slp;
    }
  for(i=0;i<idata;i++){
     itcpt=datay[i]-datax[i]*slpmx;
     if(itcptmn>itcpt)itcptmn=itcpt;
     itcpt=datay[i]-datax[i]*slpmn;
     if(itcptmx<itcpt)itcptmx=itcpt;
    }

 // nstep=atoi(argv[3]);

  for(step=0;;step++){
//    printf("slpmin: %f  slpmax: %f\n",slpmn,slpmx);
    slpstp=(slpmx-slpmn)/9;
    itcptstp=(itcptmx-itcptmn)/9;

    rsdmn=999999999.;
    for(i=0;i<10;i++){
       slp=slpmn+i*slpstp;
       for(j=0;j<10;j++){
          itcpt=itcptmn+j*itcptstp;
          rsd[i*10+j]=0;
          slpp[i*10+j]=slp;
          itcptt[i*10+j]=itcpt;
          for(ii=0;ii<idata;ii++){
             rsd[i*10+j]+=(slp*datax[ii]+itcpt-datay[ii])*(slp*datax[ii]+itcpt-datay[ii]);
            }
          rsd[i*10+j]=sqrt(rsd[i*10+j]/(idata-1));
         }
      }
    for(i=0;i<5;i++)
       for(j=i;j<100;j++)
          if(rsd[i]>rsd[j]){
             rsdfit=rsd[i];
             slpfit=slpp[i];
             itcptfit=itcptt[i];
             rsd[i]=rsd[j];
             slpp[i]=slpp[j];
             itcptt[i]=itcptt[j];
             rsd[j]=rsdfit;
             slpp[j]=slpfit;
             itcptt[j]=itcptfit;
            }

    slpmn=999999999;slpmx=-999999999;
    itcptmn=999999999;itcptmx=-999999999;
    for(i=0;i<5;i++){
//printf("residue: %f\n",rsd[i]);
       if(slpmn>slpp[i])slpmn=slpp[i];
       if(slpmx<slpp[i])slpmx=slpp[i];
       if(itcptmn>itcptt[i])itcptmn=itcptt[i];
       if(itcptmx<itcptt[i])itcptmx=itcptt[i];
      }
    if((slpmn-slpmx)*(slpmn-slpmx)<slpmn*slpmn/100000)break;
   }
//  printf("steps: %d\n",step);
//  printf("slope: %f  intercept: %f  residue: %f\n",slpp[0],itcptt[0],rsd[0]);
//  printf("%f %f %f\n",slpp[0],itcptt[0],rsd[0]);
  *slpout=slpp[0]; *itcptout=itcptt[0]; *rsdout=rsd[0];
  return 1;
}


int gau_fit_line ( double *datax, double *datay, int idata, double slpmn, double slpmx, double xmid, double alpha, double *slpout, double *itcptout, double *rsdout )
{
  FILE *f;
  int i,j,ii;
  double weit, weight;
  double slp,slpstp,itcpt,itcptmx,itcptmn,itcptstp,slpfit,itcptfit,slpp[100],itcptt[100],rsd[100],dataxx,datayy,datamin,rsdmn,rsdfit;

  itcptmx=-999999999.;itcptmn=999999999.;
  for(i=0;i<idata;i++){
     itcpt=datay[i]-datax[i]*slpmx;
     if(itcptmn>itcpt)itcptmn=itcpt;
     itcpt=datay[i]-datax[i]*slpmn;
     if(itcptmx<itcpt)itcptmx=itcpt;
    }


  for(;;){
//    printf("slpmin: %f  slpmax: %f\n",slpmn,slpmx);
    slpstp=(slpmx-slpmn)/9;
    itcptstp=(itcptmx-itcptmn)/9;

    rsdmn=999999999.;
    for(i=0;i<10;i++){
       slp=slpmn+i*slpstp;
       for(j=0;j<10;j++){
          itcpt=itcptmn+j*itcptstp;
          rsd[i*10+j]=0; weight=0;
          slpp[i*10+j]=slp;
          itcptt[i*10+j]=itcpt;
          for(ii=0;ii<idata;ii++){
             weit=exp(-alpha*pow((datax[ii]-xmid),2)); weight+=weit;
             rsd[i*10+j] += pow((datay[ii]-itcpt-datax[ii]*slp),2)*weit;
             //rsd[i*10+j]+=(slp*datax[ii]+itcpt-datay[ii])*(slp*datax[ii]+itcpt-datay[ii]);
            }
          //rsd[i*10+j]=sqrt(rsd[i*10+j]/(idata-1));
          rsd[i*10+j] = sqrt(rsd[i*10+j])/weight;
         }
      }
    for(i=0;i<5;i++)
       for(j=i;j<100;j++)
          if(rsd[i]>rsd[j]){
             rsdfit=rsd[i];
             slpfit=slpp[i];
             itcptfit=itcptt[i];
             rsd[i]=rsd[j];
             slpp[i]=slpp[j];
             itcptt[i]=itcptt[j];
             rsd[j]=rsdfit;
             slpp[j]=slpfit;
             itcptt[j]=itcptfit;
            }

    slpmn=999999999;slpmx=-999999999;
    itcptmn=999999999;itcptmx=-999999999;
    for(i=0;i<5;i++){
//printf("residue: %f\n",rsd[i]);
       if(slpmn>slpp[i])slpmn=slpp[i];
       if(slpmx<slpp[i])slpmx=slpp[i];
       if(itcptmn>itcptt[i])itcptmn=itcptt[i];
       if(itcptmx<itcptt[i])itcptmx=itcptt[i];
      }
    if((slpmn-slpmx)*(slpmn-slpmx)<slpmn*slpmn/100000)break;
   }
//  printf("slope: %f  intercept: %f  residue: %f\n",slpp[0],itcptt[0],rsd[0]);
//  printf("%f %f %f\n",slpp[0],itcptt[0],rsd[0]);
  *slpout=slpp[0]; *itcptout=itcptt[0]; *rsdout=rsd[0];
  return 1;
}

int gaus_smooth ( double *datax, double *datay, int idata, double slpmn, double slpmx, double dx, double *smthy, double *rsdout, int *ism )
{
  int i, j, ii, ib, ie, nstep=10;
  double dataxx, datayy, datamin, hwidth=100., alpha=0.5/hwidth/hwidth, rsdmx=0.06;
  double slp, itcpt, step, rsdmin;

  for(i=0;i<idata;i++) {
     datamin=datax[i]; ii=i;
     for(j=i+1;j<idata;j++)
        if(datamin>datax[j]){
           datamin=datax[j];
           ii=j;
        }
     if(ii==i) continue;
     dataxx=datax[i];
     datayy=datay[i];
     datax[i]=datax[ii];
     datay[i]=datay[ii];
     datax[ii]=dataxx;
     datay[ii]=datayy;
  }

  double xtmp, slpmin, slpmax, slpb, itcptb;
  int iseg, segb, sege, dxstep=(int)ceil(datax[idata-1]/dx)+2;
  *ism = dxstep;
  sege=-1; slpb=(slpmn+slpmx)/2.;
  for(iseg=0;;) {
     if(sege==dxstep) break;
     slp = slpb;
     rsdmin=99999.;
     segb=-1;
     for(i=sege+1;i<dxstep;i++) {
        xtmp=i*dx;
        if(xtmp<datax[0]-2.*dx) continue;
        for(j=idata-1;j>=0;j--) if(datax[j]<xtmp-2*hwidth) break;
        ib=j+1;
        for(j=0;j<idata;j++) if(datax[j]>xtmp+2*hwidth) break;
        ie=j-1;
        if(ib > ie)  { smthy[i] = -99999.; rsdout[i]=-999.; }
        else if(ib == ie) { smthy[i] = datay[ib]+slp*(xtmp-datax[ib]); rsdout[i] = -rsdmx*3.; }
        else if(ib == ie-1) {smthy[i] = datay[ib]+(datay[ie]-datay[ib])*(xtmp-datax[ib])/dx; rsdout[i]=-rsdmx*2.;}
        else gau_fit_line ( &datax[ib], &datay[ib], ie-ib+1, slpmn, slpmx, xtmp, alpha, &slp, &itcpt, &rsdout[i] );
        if(rsdout[i]>rsdmx) { smthy[i] = slp*xtmp+itcpt; rsdout[i]=-rsdmx*2.; }
        if(segb==-1) { if(rsdout[i]>=0) segb=i; else continue; }
        else if(rsdout[i]<0) { sege=i; break; }
        //cout<<"rsd: "<<rsdout[i]<<endl;
        if(rsdout[i]<rsdmin) {rsdmin=rsdout[i]; ii=i; slpb=slp; itcptb=itcpt;}
        //else rsdout[i] /= sqrt(ie-ib+1);
     }
     if(i==dxstep) {
        if(segb==-1) break;
        else sege=i;
     }
     //cout<<iseg<<" "<<segb*dx<<" to "<<sege*dx<<endl;
     //cout<<"best at: "<<ii*dx<<"  slpb: "<<slpb<<endl;
     slp = slpb;
     smthy[ii] = slpb*ii*dx+itcptb;
     for(i=ii+1;i<sege;i++) {
        xtmp=i*dx;
        if(xtmp<datax[0]-2.*dx) continue;
        for(j=idata-1;j>=0;j--) if(datax[j]<xtmp-2*hwidth) break;
        ib=j+1;
        for(j=0;j<idata;j++) if(datax[j]>xtmp+2*hwidth) break;
        ie=j-1;
        if(ib > ie-2) continue;
        slpmax = (slpmx-slpmn)/(100./dx);
        slpmin=slp-slpmax ; slpmax=slp+slpmax;
        if(slpmin<slpmn) slpmin=slpmn;
        if(slpmax>slpmx) slpmax=slpmx;
        gau_fit_line ( &datax[ib], &datay[ib], ie-ib+1, slpmin, slpmax, xtmp, alpha, &slp, &itcpt, &rsdout[i] );
      //cout<<xtmp<<" "<<slp<<endl;
        smthy[i] = slp*xtmp+itcpt;
     }
     slp = slpb;
     for(i=ii-1;i>=segb;i--) {
        xtmp=i*dx;
        if(xtmp<datax[0]-2.*dx) continue;
        for(j=idata-1;j>=0;j--) if(datax[j]<xtmp-2*hwidth) break;
        ib=j+1;
        for(j=0;j<idata;j++) if(datax[j]>xtmp+2*hwidth) break;
        ie=j-1;
        if(ib > ie-2) continue;
        slpmax = (slpmx-slpmn)/(100./dx);
        slpmin=slp-slpmax ; slpmax=slp+slpmax;
        if(slpmin<slpmn) slpmin=slpmn;
        if(slpmax>slpmx) slpmax=slpmx;
        gau_fit_line ( &datax[ib], &datay[ib], ie-ib+1, slpmin, slpmax, xtmp, alpha, &slp, &itcpt, &rsdout[i] );
      //cout<<xtmp<<" "<<slp<<endl;
        smthy[i] = slp*xtmp+itcpt;
     }
     if((sege-segb)*dx<1.5*hwidth) { for(i=segb;i<sege;i++) rsdout[i]=-rsdmx*2.; continue; }
     iseg++;
  }
  cout<<iseg<<" segments smoothed"<<endl;

  return 1;
}


#define NSTA 2000
int main (int argc, char *argv[])
{
if (argc != 6) {
  cout<<"please input [amp_file] [beta_file] [pvel_file] [station.lst] [per]"<<endl;
  return 0;
  }
   int i,j,ii,nsta,nlst;
   int ilon, ilat;
   char buff[300];
   FILE *inf, *ouf;
   float DGR = 0.1;
   double lonmin=360, lonmax=0, latmin=90, latmax=-90, per=atof(argv[5]);
   double lon[NSTA], lat[NSTA], azi[NSTA], dis[NSTA], amp[NSTA], clong, clati, lontmp, lattmp;
   char stalst[NSTA][10];
   double lonlst[NSTA], latlst[NSTA];
   double dattmp;

   //clong=atof(argv[4]); clati=atof(argv[5]);
   if(clong<0) clong += 360.;

   if((inf = fopen(argv[4],"r")) == NULL) {
      cout<<"Can't open station list "<<argv[6]<<endl;
      return 0;
   }
   for(i=0;;i++){
      if((fgets(buff, 300, inf)) == NULL)break;
      sscanf(buff,"%s %lf %lf", &stalst[i][0], &lonlst[i], &latlst[i]);
      if(lonlst[i]<0) lonlst[i] += 360.;
   }
   nlst=i;
   fclose(inf);

   if((inf = fopen(argv[1],"r")) == NULL) {
      cout<<"Can't open amp file "<<argv[1]<<endl;
      return 0;
   }
   fgets(buff, 300, inf);
   sscanf(buff,"%lf %lf", &clong, &clati);
   for(i=0;;i++){
      if((fgets(buff, 300, inf)) == NULL)break;
      sscanf(buff,"%lf %lf %lf", &lon[i], &lat[i], &amp[i]);
      if(lon[i]<0) lon[i] += 360.;
      calc_azimuth(clati,clong,lat[i],lon[i],&azi[i]);
      dis[i]=get_dist(clati, clong, lat[i], lon[i]);
      if(lonmin>lon[i]) lonmin = lon[i];
      if(lonmax<lon[i]) lonmax = lon[i];
      if(latmin>lat[i]) latmin = lat[i];
      if(latmax<lat[i]) latmax = lat[i];
   }
   nsta=i;
   fclose(inf);

   lonmin = floor(lonmin/DGR)*DGR; lonmax = ceil(lonmax/DGR)*DGR;
   latmin = floor(latmin/DGR)*DGR; latmax = ceil(latmax/DGR)*DGR;
   int nptslon=(int)((lonmax-lonmin)/DGR)+1;
   int nptslat=(int)((latmax-latmin)/DGR)+1;
   double beta[nptslon][nptslat], pvel[nptslon][nptslat], cmean; //ydata[nptslon][nptslat], azi[nptslon][nptslat], dist[nptslon][nptslat];

   for(i=0;i<nptslon;i++) for(j=0;j<nptslat;j++) {
      beta[i][j] = -1; pvel[i][j] = -1; }
/*
   for(i=0;i<nptslon;i++) for(j=0;j<nptslat;j++) {
      lontmp = lonmin+i*DGR; lattmp = latmin+j*DGR;
      dist[i][j] = get_dist(clati, clong, lattmp, lontmp);
      calc_azimuth(clati,clong,lattmp,lontmp,&azi[i][j]);
   }
*/

   ouf = fopen("region_tmp", "w");
   fprintf(ouf, "-R%.0f/%.0f/%.0f/%.0f\n", floor(lonmin), ceil(lonmax), floor(latmin), ceil(latmax));
   fclose(ouf);

/*
   sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel_T0.2 %s region_tmp %f\n", argv[1], DGR);
   system(buff);
   sprintf(buff,"%s.HD_0.2\0", argv[1]);
   if((inf = fopen(buff,"r")) == NULL) {
      cout<<"Can't open amp HD file "<<buff<<endl;
      return 0;
   }
   for(i=0;;i++){
      if((fgets(buff, 300, inf)) == NULL)break;
      sscanf(buff,"%lf %lf %lf", &lontmp, &lattmp, &dattmp);
      ilon=(int)floor((lontmp-lonmin)/DGR+0.5);
      ilat=(int)floor((lattmp-latmin)/DGR+0.5);
      if(ilon<0 || ilon>=nptslon || ilat<0 || ilat>=nptslat) continue;
      ampg[ilon][ilat] = dattmp;
   }
   fclose(inf);
*/

   sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel %s region_tmp %f\n", argv[2], DGR);
   system(buff);
   sprintf(buff,"%s.HD\0", argv[2]);
   if((inf = fopen(buff, "r")) == NULL) {
      cout<<"Can't open beta HD file "<<buff<<endl;
      return 0;
   }
   for(i=0;;i++){
      if((fgets(buff, 300, inf)) == NULL)break;
      sscanf(buff,"%lf %lf %lf", &lontmp, &lattmp, &dattmp);
      ilon=(int)floor((lontmp-lonmin)/DGR+0.5);
      ilat=(int)floor((lattmp-latmin)/DGR+0.5);
      if(ilon<0 || ilon>=nptslon || ilat<0 || ilat>=nptslat) continue;
      beta[ilon][ilat] = dattmp;
   }
   fclose(inf);

   sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel %s region_tmp %f\n", argv[3], DGR);
   system(buff);
   sprintf(buff,"%s.HD\0", argv[3]);
   if((inf = fopen(buff,"r")) == NULL) {
      cout<<"Can't open pvel HD file "<<buff<<endl;
      return 0;
   }
   cmean=0; j=0;
   for(i=0;;i++){
      if((fgets(buff, 300, inf)) == NULL)break;
      sscanf(buff,"%lf %lf %lf", &lontmp, &lattmp, &dattmp);
      ilon=(int)floor((lontmp-lonmin)/DGR+0.5);
      ilat=(int)floor((lattmp-latmin)/DGR+0.5);
      if(ilon<0 || ilon>=nptslon || ilat<0 || ilat>=nptslat) continue;
      pvel[ilon][ilat] = dattmp;
      cmean += dattmp;
      j++;
   }
   fclose(inf);
   cmean /= j;

/*
   for(i=0;i<nptslon;i++) for(j=0;j<nptslat;j++) {
      if(ampg[i][j] == -1 || beta[i][j] == -1 || pvel[i][j] == -1) {
         ydata[i][j]=-1; continue; }
//      cout<<pvel[i][j]<<" "<<beta[i][j]<<" "<<log(1./beta[i][j]*100.*sqrt(1./pvel[i][j]*3.5))<<endl;
      ydata[i][j] = log(ampg[i][j]/beta[i][j]*100.*sqrt(dist[i][j]/pvel[i][j]*3.5));
      ampg[i][j] = log(ampg[i][j]*sqrt(dist[i][j]));
   }
*/

   int bin=5, nbin=(int)(359.9/bin)+1;
   int ibin, iazi, idat, ndat[nbin];
   short datn[nbin][300], datan[nbin][300];
   char datsta[nbin][300][10], datasta[nbin][500][10];
   double daty[nbin][300], datx[nbin][300], datay[nbin][500], datax[nbin][500], smthy[nbin][500], smrsd[500];
   double slp, itcpt[nbin], itcpt2[nbin], ditcpt, rsd[nbin], tazi, tdis, dismin, dismax;
   double pi=3.1415926536, slpmax=-2.*pi/per/cmean/600., slpmin=-2.*pi/per/cmean/20.;

   for(ibin=0;ibin<nbin;ibin++){
      iazi=ibin*bin;
      idat=0;
      for(i=0;i<nsta;i++) {
         tazi=azi[i];
         if( (tazi >= iazi+bin && tazi < iazi-bin+360.) || (tazi < iazi-bin && tazi >= iazi+bin-360.)) continue;
         for(j=0;j<nlst;j++) 
            if((fabs(lat[i]-latlst[j])+fabs(lon[i]-lonlst[j]))<0.01) break;
         if(j==nlst) sprintf(datsta[ibin][idat],"NaN\0");
         else sprintf(datsta[ibin][idat],"%s\0",stalst[j]);
         tdis=dis[i];
         ilon=(int)floor((lon[i]-lonmin)/DGR+0.5);
         ilat=(int)floor((lat[i]-latmin)/DGR+0.5);
         //datyold[idat] = log(amp[i]*sqrt(tdis));
         daty[ibin][idat] = log(amp[i]/beta[ilon][ilat]*100.*sqrt(tdis/pvel[ilon][ilat]*3.5));
         datx[ibin][idat] = tdis;
         datn[ibin][idat] = i;
         idat++;
//         if(dismin>tdis) dismin=tdis;
//         if(dismax<tdis) dismax=tdis;
      }
      if(idat<6) { itcpt[ibin]=-999; continue; }
      fit_line ( &datx[ibin][0], &daty[ibin][0], idat, &slp, &itcpt[ibin], &rsd[ibin] );
      //cout<<"azi: "<<iazi<<"  nsta: "<<idat<<"  slp: "<<slp*1e3<<"  itcpt: "<<itcpt[ibin]<<"  rsd: "<<rsd[ibin]<<endl;
      for(i=0;i<idat;i++) {
         if(fabs(daty[ibin][i]-itcpt[ibin]-slp*datx[ibin][i])>rsd[ibin]*1.5) {
           for(j=i+1;j<idat;j++) {
              sprintf(datsta[ibin][j-1],datsta[ibin][j]);
              datx[ibin][j-1]=datx[ibin][j]; daty[ibin][j-1]=daty[ibin][j];
              datn[ibin][j-1]=datn[ibin][j];
           }
           idat--;
         }
      }
      if(idat<6) { itcpt[ibin]=-999; continue; }
      fit_line ( &datx[ibin][0], &daty[ibin][0], idat, &slp, &itcpt[ibin], &rsd[ibin] );
      if(slp>slpmax || slp<slpmin || rsd[ibin]>0.6) { itcpt[ibin]=-999; continue; }
      ndat[ibin] = idat;
/*
      sprintf(buff,"alpha_r_curve_fit_%ddeg\0", iazi);
      ouf=fopen(buff,"w");
      for(i=0;i<idat;i++) fprintf(ouf, "%8.3f\t%7.4f\t%7.4f\t%7.4f\t%s\n", datx[ibin][i], daty[ibin][i], itcpt[ibin]+slp*datx[ibin][i], rsd[ibin]*1.5, datsta[ibin][i]);
      fclose(ouf);
*/
   }

   int ibinl, ibinh, ism;
   double ddis=10.;
   double grdtl, grdth, alpha[nsta], weight[nsta];
   for(i=0;i<nsta;i++) { alpha[i]=0.; weight[i]=0.; }
   for(ibin=0;ibin<nbin;ibin++){
      ibinl=ibin-1; ibinh=ibin+1;
      if(ibinl<0) ibinl += nbin;
      if(ibinh>=nbin) ibinh -= nbin;
      if(itcpt[ibinh] == -999 || itcpt[ibinl] == -999) continue;
      iazi=ibin*bin;
      ditcpt = itcpt[ibinh] - itcpt[ibinl];
      for(i=0;i<ndat[ibinl];i++) {
         sprintf(datasta[ibin][i],"%s\0",datsta[ibinl][i]);
         datax[ibin][i]=datx[ibinl][i]; datay[ibin][i]=daty[ibinl][i]+ditcpt;
         datan[ibin][i]=datn[ibinl][i];
      }
      for(j=0;j<ndat[ibinh];j++) {
         sprintf(datasta[ibin][i+j],"%s\0",datsta[ibinh][j]);
         datax[ibin][i+j]=datx[ibinh][j]; datay[ibin][i+j]=daty[ibinh][j];
         datan[ibin][i+j]=datn[ibinh][j];
      }
      idat=i+j;

      fit_line ( &datax[ibin][0], &datay[ibin][0], idat, &slp, &itcpt2[ibin], &rsd[ibin] );
      for(i=0;i<idat;i++) {
         if(fabs(datay[ibin][i]-itcpt2[ibin]-slp*datax[ibin][i])>rsd[ibin]*1.5) {
           for(j=i+1;j<idat;j++) {
              sprintf(datasta[ibin][j-1],datasta[ibin][j]);
              datax[ibin][j-1]=datax[ibin][j]; datay[ibin][j-1]=datay[ibin][j];
           }
           idat--;
         }
      }
  //if(iazi<120) continue;
      cout<<iazi<<endl;
      gaus_smooth ( &datax[ibin][0], &datay[ibin][0], idat, slpmin, slpmax, ddis, &smthy[ibin][0], &smrsd[0], &ism ); 
  //exit(0);

      for(i=0;i<idat;i++) {
         ii= int(datax[ibin][i]/ddis);
         //if(smrsd[ii-1]>0.06 || smrsd[ii]>0.06 || smrsd[ii+1]>0.06 || smrsd[ii+2]>0.06 ) continue;
         if(smrsd[ii-1]<0 || smrsd[ii]<0 || smrsd[ii+1]<0 || smrsd[ii+2]<0) continue;
         grdtl = (smthy[ibin][ii+1]-smthy[ibin][ii-1])/2./ddis;
         grdth = (smthy[ibin][ii+2]-smthy[ibin][ii])/2./ddis;
         j=datan[ibin][i];
         alpha[j] += -(grdtl+(grdth-grdtl)*(datax[ibin][i]-ii*ddis)/ddis)/fabs(iazi-azi[j]+1e-10);
         weight[j] += 1./fabs(iazi-azi[j]+1e-10);
      }

      sprintf(buff,"alpha_r_sta_%ddeg\0", iazi);
      ouf=fopen(buff,"w");
      for(i=0;i<idat;i++) fprintf(ouf, "%8.3f\t%7.4f\t%s\n", datax[ibin][i], datay[ibin][i], datasta[ibin][i]);
      fclose(ouf);
      sprintf(buff,"alpha_r_curve_fit_%ddeg\0", iazi);
      ouf=fopen(buff,"w");
      for(i=0;i<ism;i++) {
         if(i*ddis<datax[ibin][0]-2.*ddis) continue;
         if(smthy[ibin][i] < -99998.) continue;
         fprintf(ouf, "%8.3f\t%7.4f\t%7.4f\n", i*ddis, smthy[ibin][i], smrsd[i]);
      }
      fclose(ouf);
    //exit(0);
   }

   sprintf(buff,"%s_Q_map\0", argv[1]);
   ouf=fopen(buff,"w");
   for(i=0;i<nsta;i++) {
      if(weight[i] < 0.2) continue;
      alpha[i] /= weight[i];
      if(alpha[i] < -slpmax || alpha[i] > -slpmin) continue;
      ilon=(int)floor((lon[i]-lonmin)/DGR+0.5);
      ilat=(int)floor((lat[i]-latmin)/DGR+0.5);
      fprintf(ouf,"%9.4f %8.4f %6.2f\n",lon[i],lat[i],2*pi/per/pvel[ilon][ilat]/alpha[i]);
   }
   fclose(ouf);
//slpmax=-2.*pi/per/cmean/600., slpmin

   return 1;
}
