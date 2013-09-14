#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

double e_n=2.7182818284590451;


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
  double slp,slpmx,slpmn,slpstp,itcpt,itcptmx,itcptmn,itcptstp,slpfit,itcptfit,slpp[100],itcptt[100],rsd[100],dataxx,datayy,rsdmn,rsdfit;


  for(i=0;i<idata;i++)
     for(j=i;j<idata;j++)
        if(datax[i]>datax[j]){
           dataxx=datax[i];
           datayy=datay[i];
           datax[i]=datax[j];
           datay[i]=datay[j];
           datax[j]=dataxx;
           datay[j]=datayy;
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
cout<<slpp[0]<<endl;
//  printf("steps: %d\n",step);
//  printf("slope: %f  intercept: %f  residue: %f\n",slpp[0],itcptt[0],rsd[0]);
//  printf("%f %f %f\n",slpp[0],itcptt[0],rsd[0]);
  *slpout=slpp[0]; *itcptout=itcptt[0]; *rsdout=rsd[0];
  return 1;
}


int  get_ave(double *input, int cvn,double *ave, double *std ) {
  int i,j,k;
  double tave = 0.;
  double tstd = 0.;
  for (i=0;i<cvn;i++) {
    tave= tave+ input[i];
    }
  tave = tave/cvn;
  for (i=0;i<cvn;i++) {
    tstd = tstd + (tave-input[i])*(tave-input[i]);
    }
  if (cvn>2) {
    tstd = sqrt(tstd/(cvn-1));
    }
  else
    tstd = 10000.;
  *ave=tave;
  *std = tstd;

  return 0;
  }

int main (int argc, char *argv[])
{
int NSTA = 100000;
int i,j,k,ii;
char outname[300], outname2[300], buff[300];
double lonb[NSTA],latb[NSTA],lonc[NSTA],latc[NSTA],beta[NSTA],pvel[NSTA];
double distb[NSTA],azib[NSTA],distc[NSTA],azic[NSTA];
double pi=3.14159265359, lon1 = -114.8757+360,lat1 = 41.1003, tazi,tdis;
double beta0=-1, pvel0=-1;
FILE *inf,*outf,*outf2,*outa1,*outa2,*outa3;
int nstab,nstac;
if (argc != 4) {
  cout<<"please input [beta_infile] [phase_infile] [azi_bin]"<<endl;
  return 0;
  }

float lonminb=360, lonmaxb=-360, latminb=90, latmaxb=-90;
float lonminc=360, lonmaxc=-360, latminc=90, latmaxc=-90;
double lontmp=int(lon1*10+0.5)/10., lattmp=int(lat1*10+0.5)/10.;
cout<<lontmp<<" "<<lattmp<<endl;
inf = fopen(argv[1],"r");
for (i=0;;i++) {
  if(fgets(buff, 300, inf) == NULL )break;
  sscanf(buff,"%lf %lf %lf", &lonb[i],&latb[i],&beta[i]);
  //beta0=beta[i];
  if((lontmp == lonb[i]) && (lattmp == latb[i])) beta0=beta[i];
  if(lonminb>lonb[i])lonminb=lonb[i];
  if(lonmaxb<lonb[i])lonmaxb=lonb[i];
  if(latminb>latb[i])latminb=latb[i];
  if(latmaxb<latb[i])latmaxb=latb[i];
  }
nstab = i;
fclose(inf);
inf = fopen(argv[2],"r");
for (i=0;;i++) {
  if (fscanf(inf,"%lf %lf %lf", &(lonc[i]),&(latc[i]),&(pvel[i])) != 3) break;
  if (lontmp==lonc[i] && lattmp==latc[i]) pvel0=pvel[i];
  if(lonminc>lonc[i])lonminc=lonc[i];
  if(lonmaxc<lonc[i])lonmaxc=lonc[i];
  if(latminc>latc[i])latminc=latc[i];
  if(latmaxc<latc[i])latmaxc=latc[i];
  }
nstac = i;
fclose(inf);
if(lonminb<lonminc) lonminb=lonminc;
if(lonmaxb>lonmaxc) lonmaxb=lonmaxc;
if(latminb<latminc) latminb=latminc;
if(latmaxb>latmaxc) latmaxb=latmaxc;

float dgr=0.1;
int lonnpts=int((lonmaxb-lonminb)/dgr)+1, latnpts=int((latmaxb-latminb)/dgr)+1;
double alpha[lonnpts][latnpts],alpha1[lonnpts][latnpts], alpha2[lonnpts][latnpts], dist[lonnpts][latnpts], azi[lonnpts][latnpts];
double glon,glat;

sprintf(outname,"alpha_beta_map\0");
outa1=fopen(outname,"w");
sprintf(outname,"alpha_pvel_map\0");
outa2=fopen(outname,"w");
sprintf(outname,"alpha_map\0");
outa3=fopen(outname,"w");
for(i=0;i<lonnpts;i++){
   for(j=0;j<latnpts;j++){
      glon=lonminb+i*dgr;
      glat=latminb+j*dgr;
      dist[i][j] = get_dist(glat,glon,lat1,lon1);
      calc_azimuth(lat1,lon1,glat,glon,&azi[i][j]);
      for(k=0;k<nstab;k++)
         if(pow(glon-lonb[k],2)<0.01 && pow(glat-latb[k],2)<0.01) {
            alpha1[i][j] = -1.e-3/dist[i][j]*log(beta[k]/beta0);
            //alpha1[i][j] = -log(beta[k]/beta0);
            break;
         }
      if(k==nstab) {alpha[i][j]=-1; continue;}
      for(k=0;k<nstac;k++) 
         if(pow(glon-lonc[k],2)<0.01 && pow(glat-latc[k],2)<0.01) {
            alpha2[i][j] = -0.5e-3/dist[i][j]*log(pvel[k]/pvel0);
            //alpha2[i][j] = -0.5*log(pvel[k]/pvel0);
            break;
         }
      if(k==nstac) {alpha[i][j]=-1; continue;}
      alpha[i][j] = alpha1[i][j] + alpha2[i][j];
      fprintf(outa3,"%.1f\t%.1f\t%g\n",glon,glat,alpha[i][j]);
      fprintf(outa1,"%.1f\t%.1f\t%g\n",glon,glat,alpha1[i][j]);
      fprintf(outa2,"%.1f\t%.1f\t%g\n",glon,glat,alpha2[i][j]);
   }
}
fclose(outa1);
fclose(outa2);
fclose(outa3);
return 2;
float azibin,azimax,azimin,step;
double talpha[NSTA],tdist[NSTA],dtmp;
int nazi,ngrd;
azibin = atof(argv[3]);
nazi = 100;
step = 360./nazi;

for (k=0;k<nazi;k++) {
  azimax = k*step + azibin/2.;
  azimin = azimax - azibin;
  ngrd = 0;
  for(i=0;i<lonnpts;i++){
     for(j=0;j<latnpts;j++){
        if ( (azi[i][j] > azimax && azi[i][j] < azimin+360.) || (azi[i][j] < azimin && azi[i][j] > azimax -360.)) continue;
        if(alpha[i][j]==-1) continue;
        talpha[ngrd] = alpha[i][j];
        tdist[ngrd] = dist[i][j];
        ngrd++;
     }//j
  }//i
  for(i=0;i<ngrd;i++)
     for(j=i;j<ngrd;j++)
        if(tdist[i]>tdist[j]){
           dtmp=tdist[i];
           tdist[i]=tdist[j];
           tdist[j]=dtmp;
           dtmp=talpha[i];
           talpha[i]=talpha[j];
           talpha[j]=dtmp;
           }
  sprintf(outname2,"%.1f_dist_alpha_model\0",k*step);
  outf2 = fopen(outname2,"w");
  for(ii=0;ii<ngrd;ii++) fprintf(outf2,"%f\t%f\n",tdist[ii],talpha[ii]);
  fclose(outf2);
}//k


/*
  fit_line( &tdist[0], &tamp[0], k, &slp, &itcpt, &std1 );
  sprintf(outname2,"%.1f_dist_amp\0",i*step);
  outf2 = fopen(outname2,"w");
  for(ii=0;ii<k;ii++) fprintf(outf2,"%f\t%f\n",tdist[ii],tamp[ii]);
  fclose(outf2);
  sprintf(outname2,"%.1f_dist_ampcrctd\0",i*step);
  outf2 = fopen(outname2,"w");
  for(ii=0;ii<k;ii++) fprintf(outf2,"%f\t%f\n",tdist[ii],exp(tamp[ii]-slp*tdist[ii]-itcpt));
  fclose(outf2);
  sprintf(outname2,"%.1f_dist_alpha\0",i*step);
  outf2 = fopen(outname2,"w");
  for(ii=0;ii<k;ii++) fprintf(outf2,"%f\t%f\n",tdist[ii],-(tamp[ii]-itcpt)/tdist[ii]);
  fclose(outf2);


//  printf("%f %f %f\n",slp,itcpt,std1);
//  ave1 = 0.;
//  std1 = 0.;
//  get_ave(&(tamp[0]),cvn,&ave1,&std1);
  if (k>2)
     fprintf(outf,"%g %g %g %g\n", i*step,exp(itcpt),std1*exp(itcpt),-slp);
//  abort();
  }
fclose(outf);

double A0,A2,phi2,alpha;
  sprintf(buff,"python /home/tianye/code/Programs/ASN/do_fourier.py %s | awk '{print $1,$3,$5}' > temp.fit",outname);
  system(buff);
  outf = fopen("temp.fit","r");
  fscanf(outf,"%lf %lf %lf", &A0, &A2,&phi2);
  fclose(outf);

  sscanf(argv[1],"%[^_]",buff);
  strcat(buff,"_am.txt_v2");
  inf = fopen(buff,"r");
  for (i=0;;i++) {
     if (fscanf(inf,"%lf %lf %lf", &(lon[i]),&(lat[i]),&(amp[i])) != 3) break;
  }
  nsta0 = i;


  strcat(buff,"_azi_norm");
  outf = fopen(buff,"w");
  for(i=0;i<nsta0;i++) {
     calc_azimuth(lat1, lon1, lat[i], lon[i], &alpha);
     amp[i] /= A0+A2*sin(pi/180.*(2*alpha+phi2));
     fprintf(outf,"%lf\t%lf\t%lf\n", lon[i],lat[i],amp[i]);
  }
  fclose(outf);
*/
  return 1; 
}
