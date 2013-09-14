#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

int gau_fit_line ( double *datax, double *datay, double *oweit, int idata, double slpmn, double slpmx, double xmid, double alpha, double *slpout, double *itcptout, double *rsdout )
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
             weit=exp(-alpha*pow((datax[ii]-xmid),2))*oweit[ii]; weight+=weit;
             rsd[i*10+j] += pow((datay[ii]-itcpt-datax[ii]*slp),2)*weit;
             //rsd[i*10+j]+=(slp*datax[ii]+itcpt-datay[ii])*(slp*datax[ii]+itcpt-datay[ii]);
            }
          //rsd[i*10+j]=sqrt(rsd[i*10+j]/(idata-1));
          rsd[i*10+j] = sqrt(rsd[i*10+j])/(weight-1);
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

int gaus_smooth ( double *datax, double *datay, double *oweit, int idata, double slpmn, double slpmx, double hwidth, double dx, double *cfdt, double *smthy, double *rsdout, int *ism )
{
  int i, j, ii, ipt, ib, ie, nstep=10;
  double dataxx, datayy, weitt, weight, datamin, alpha=0.5/hwidth/hwidth, rsdmx=1e3;
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
     weitt=oweit[i];
     datax[i]=datax[ii];
     datay[i]=datay[ii];
     oweit[i]=oweit[ii];
     datax[ii]=dataxx;
     datay[ii]=datayy;
     oweit[ii]=weitt;
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
        if(ib>ie || datax[ib]>xtmp+hwidth || datax[ie]<xtmp-hwidth) { smthy[i] = 0.; rsdout[i]=-999.; }
        //if(ib > ie)  { smthy[i] = 0.; rsdout[i]=-999.; }
        else if(ib == ie) { smthy[i] = datay[ib]+slp*(xtmp-datax[ib]); rsdout[i] = -rsdmx*3.; }
        else if(ib == ie-1) {smthy[i] = datay[ib]+(datay[ie]-datay[ib])*(xtmp-datax[ib])/dx; rsdout[i]=-rsdmx*2.;}
        else gau_fit_line ( &datax[ib], &datay[ib], &oweit[ib], ie-ib+1, slpmn, slpmx, xtmp, alpha, &slp, &itcpt, &rsdout[i] );
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
     cout<<"iseg "<<iseg<<": "<<segb*dx<<" to "<<sege*dx<<endl;
     cout<<"best at: "<<ii*dx<<"  slpb: "<<slpb<<endl;
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
        gau_fit_line ( &datax[ib], &datay[ib], &oweit[ib], ie-ib+1, slpmin, slpmax, xtmp, alpha, &slp, &itcpt, &rsdout[i] );
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
        slpmax = (slpmx-slpmn)/(hwidth*4/dx);
        slpmin=slp-slpmax ; slpmax=slp+slpmax;
        if(slpmin<slpmn) slpmin=slpmn;
        if(slpmax>slpmx) slpmax=slpmx;
        gau_fit_line ( &datax[ib], &datay[ib], &oweit[ib], ie-ib+1, slpmin, slpmax, xtmp, alpha, &slp, &itcpt, &rsdout[i] );
      //cout<<xtmp<<" "<<slp<<endl;
        smthy[i] = slp*xtmp+itcpt;
     }
     if((sege-segb)*dx<1.5*hwidth) { for(i=segb;i<sege;i++) rsdout[i]=-rsdmx*2.; continue; }
     iseg++;
  }
 
  for(i=0;i<idata;i++) {
     cfdt[i]=0; weight=0;
     for(j=i-1;j>=0;j--) if(datax[j]<datax[i]-2*hwidth) break;
     ib=j+1;
     for(j=i+1;j<idata;j++) if(datax[j]>datax[i]+2*hwidth) break;
     ie=j-1;
     ipt=0;
     for(j=ib;j<=ie;j++) {
        ii=int(datax[j]/dx+0.5);
        if(rsdout[ii]<0) continue;
        weitt = exp(-alpha*pow((datax[j]-datax[i]),2))*oweit[j];
        weight += weitt;
        cfdt[i] += pow((datay[j]-smthy[ii]),2)*weitt;
        ipt++;
     }
     if(ipt<3 || weight<1) { cfdt[i]=0; continue; }
     cfdt[i] = sqrt(cfdt[i])/(weight-1.);
     cfdt[i] = 1./cfdt[i]*oweit[i];
  }
  //cout<<iseg<<" segments smoothed"<<endl;

  return 1;
}



int main (int argc, char *argv[])
{
  if (argc != 3) {
  cout<<"Gauss_Smoothing_curve [input_file] [half_width]"<<endl;
  return 0;
  }

  FILE *ff;
  int i, idata, ism;
  char buff[300];
  double datax[5000], datay[5000], oweit[5000];
  double smthy[5000], cfdt[5000], rsdout[5000], hwidth=atof(argv[2]), dx=hwidth/2.;
  double slpmn=-1e3, slpmx=1e3;

  if((ff=fopen(argv[1],"r"))==NULL) {
     cout<<"Can't open file "<<argv[1]<<endl;
     return 0;
  }

  for(idata=0;;idata++) {
     if((fgets(buff, 300, ff))==NULL) break;
     sscanf(buff,"%lf %lf", &datax[idata], &datay[idata]);
     oweit[idata]=1;
  }
  fclose(ff);

  gaus_smooth ( datax, datay, oweit, idata, slpmn, slpmx, hwidth, dx, cfdt, smthy, rsdout, &ism );

  sprintf(buff,"%s_sm\0",argv[1]);
  ff=fopen(buff,"w");
  for(i=0;i<ism;i++) {
     if(i*dx<datax[0] || i*dx>datax[idata-1]) continue;
     fprintf(ff,"%lf %lf %lf\n", i*dx, smthy[i], rsdout[i] );
  }
//cout<<i*dx<<" "<<smthy[i]<<" "<<rsdout[i]<<" "<<cfdt[i]<<endl;
  fclose(ff);
  sprintf(buff,"%s_cfd\0",argv[1]);
  ff=fopen(buff,"w");
  for(i=0;i<idata;i++) fprintf(ff,"%lf %lf\n", datax[i], cfdt[i]);
  fclose(ff);

  return 1;
}
