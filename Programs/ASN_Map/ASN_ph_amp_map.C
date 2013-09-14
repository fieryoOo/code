#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
using namespace std;

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


int main (int argc, char *argv[])
{
  
  FILE *fin, *fin2, *fout;
  char filename[300], outname[300], buff[300], buff2[300];
  char osta[5000][10], stalst[5000][10], cstalst[5000][10], *tmp, *staname;
  char sta_pair[1000000][50], sta1[1000000][10], sta2[1000000][10];
  int i, j, k, iper, jper, nfile, nsta, ncsta, nost, nstmp;
  float tlength[1000000], ampst[5000], phst[5000], time[5000];
  float longlst[5000], latilst[5000], clonglst[5000], clatilst[5000];
  float per, per2, pertemp, perl, perh, pertmp[30], amptmp[30], phtmp[30];
  float period, amppos, snrpos, ampneg, snrneg, grvel, phvel, csamp; 
  float olong[5000], olati[5000], longtmp, latitmp, fstep=-0.08;
  double dist, alpha[5000], alphatmp, powe;
  double e_n=2.7182818284590451;
//  int flag[3600], flagl[3600]={0}, flagh[3600]={0};
  if(argc != 6) 
    {
      printf("Usage: ASN_ph_amp_map [center_sta.lst] [sta.lst] [file_time-length.lst] [period] [time-length threshold]\n"); 
      exit(-1);
    }
  
  per=atof(argv[4]);
  perl=1./exp(log(1./per)-fstep);
  perh=1./exp(log(1./per)+fstep);
  if((fin = fopen(argv[3], "r"))==NULL) {
    printf("Cannot open file list %s\n",argv[3]);
    exit(1);
  }
  for(i=0;;i++) {
    if( (fgets(buff, 100, fin)) == NULL ) break;
    fscanf(fin,"%s %f",sta_pair[i], &tlength[i]);
    tmp = strdup(sta_pair[i]);
    strtok(tmp, "_");
    staname = strtok(NULL, "_");
    strcpy(sta1[i],staname);
    staname = strtok(NULL, ".");
    strcpy(sta2[i],staname);
  }
  nfile = i;
  cout<<nfile<<" station pairs read in from the file_time-length.lst."<<endl;
  fclose(fin);

  if((fin = fopen(argv[1], "r"))==NULL) {
    printf("Cannot open center-station file %s\n",argv[1]);
    exit(1);
  }
  for(i=0;;i++) {
    if((fscanf(fin,"%s %f %f", cstalst[i], &clonglst[i], &clatilst[i]))!=3) break;
    if(clonglst[i]<0) clonglst[i]+=360;
  }
  ncsta = i;
  cout<<ncsta<<" stations read in. from center-station.lst"<<endl;
  fclose(fin);

  if((fin = fopen(argv[2], "r"))==NULL) {
    printf("Cannot open station file %s\n",argv[2]);
    exit(1);
  }
  for(i=0;;i++) {
    if((fscanf(fin,"%s %f %f", stalst[i], &longlst[i], &latilst[i]))!=3) break;
    if(longlst[i]<0) longlst[i]+=360;
  }
  nsta = i;
  cout<<nsta<<" stations read in. from station.lst"<<endl;
  fclose(fin);

  sprintf(buff,"mkdir -p Ph_Amp_Map_%ssec",argv[4]);
  system(buff);
  for(i=0;i<ncsta;i++) {
//  for(i=1298;i<1299;i++) {
//     memset(amp_n,0,3600*2*sizeof(double));
     nost=0; nstmp=0; csamp=0;
     for(j=0;j<nfile;j++) {
        if(tlength[j]<atof(argv[5])) continue;
        if(strcmp(cstalst[i],sta1[j])==0) {
           strcpy(osta[nost],sta2[j]);
           for(k=0;k<nsta;k++)
              if(strcmp(osta[nost],stalst[k])==0) {
                 olong[nost]=longlst[k]; olati[nost]=latilst[k];
                 break;
              }
           if(k==nsta) {
//             cout<<"Can't find staion "<<osta[nost]<<" in the station.lst. Skipped"<<endl;
             continue;
           }
//           calc_dist( latilst[i], longlst[i], olati[nost], olong[nost], &dist );
//           if( dist < 100 || dist > 700 ) continue;
           sprintf(filename,"%s_amp_snr\0",sta_pair[j]);
           if((fin = fopen(filename, "r"))==NULL) {
//             printf("Amp file %s not found. Skipped\n",filename);
             continue;
           }
           sprintf(filename,"%s_2_DISP.1\0",sta_pair[j]);
           if((fin2 = fopen(filename, "r"))==NULL) {
//             printf("DISP file %s not found. Skipped\n",filename);
             continue;
           }
           for(k=0;;){
              if( fgets(buff, 300, fin) == NULL ) break;
              if( fgets(buff2, 300, fin2) == NULL ){
                cout<<"amp and DISP file mismatch!"<<endl; break;}
              if((sscanf(buff,"%f %f %f %f %f", &period, &amppos, &snrpos, &ampneg, &snrneg))!=5) { 
                 cout<<"Wrong amp file format! Skipped: "<<filename<<endl; 
                 break; 
              }
              if( period < perl ) continue;
              if( period > perh ) break;
              if( snrpos>8 || (snrpos>5 && snrneg>8) ) {
                 sscanf(buff2,"%d %f %f %f %f", &iper, &pertemp, &per2, &grvel, &phvel);
                 if(per2!=period){
                   cout<<"per mismatch between amp and DISP! "<<per2<<" "<<per<<endl; break;}
                 pertmp[k] = (period-per)*(period-per)+1e-10;
                 amptmp[k] = amppos/tlength[j];
                 phtmp[k] = phvel;
                 k++;
              }
           }
           fclose(fin);
           fclose(fin2);
       
// cout<<stalst[i]<<"  "<<osta[nost]<<endl;
        }
        else if(strcmp(cstalst[i],sta2[j])==0) {
           strcpy(osta[nost],sta1[j]);
           for(k=0;k<nsta;k++)
              if(strcmp(osta[nost],stalst[k])==0) {
                 olong[nost]=longlst[k]; olati[nost]=latilst[k];
                 break;
              }
              if(k==nsta) {
//                cout<<"Can't find staion "<<osta[nost]<<" in the station.lst. Skipped"<<endl;
                continue;
              }
//           calc_dist( latilst[i], longlst[i], olati[nost], olong[nost], &dist );
//           dist_azimuth( latilst[i], longlst[i], olati[nost], olong[nost], &dist, &alpha[nost]);
//           if( dist < 100 || dist > 700 ) continue;
           sprintf(filename,"%s_amp_snr\0",sta_pair[j]);
           if((fin = fopen(filename, "r"))==NULL) {
//             printf("Amp file %s not found. Skipped\n",filename);
             continue;
           }
           sprintf(filename,"%s_2_DISP.1\0",sta_pair[j]);
           if((fin2 = fopen(filename, "r"))==NULL) {
//             printf("DISP file %s not found. Skipped\n",filename);
             continue;
           }
           for(k=0;;){
              if( fgets(buff, 300, fin) == NULL ) break;
              if( fgets(buff2, 300, fin2) == NULL ){
                cout<<"amp and DISP file mismatch!"<<endl; break;}
              if((sscanf(buff,"%f %f %f %f %f", &period, &amppos, &snrpos, &ampneg, &snrneg))!=5) { 
                 cout<<"Wrong format! Skipped: "<<filename<<endl;
                 break; }
              if( period < perl ) continue;
              if( period > perh ) break;
              if( snrneg>8 || (snrpos>8 && snrneg>5) ) {
                 sscanf(buff2,"%d %f %f %f %f", &iper, &pertmp, &per2, &grvel, &phvel);
                 if(per2!=period){
                   cout<<"per mismatch between amp and DISP!"<<per2<<" "<<per<<endl; break;}
                 pertmp[k] = (period-per)*(period-per)+1e-10;
                 amptmp[k] = ampneg/tlength[j];
                 phtmp[k] = phvel;
                 k++;
              }
           }
           fclose(fin);
           fclose(fin2);
// cout<<stalst[i]<<"  "<<osta[nost]<<endl;
        }
        else continue;

        if( k<2 ) continue;
        for(iper=0;iper<2;iper++)
           for(jper=iper+1;jper<k;jper++) {
              if(pertmp[iper]>pertmp[jper]) {
                 period=pertmp[iper];
                 pertmp[iper]=pertmp[jper];
                 pertmp[jper]=period;
                 amppos=amptmp[iper];
                 amptmp[iper]=amptmp[jper];
                 amptmp[jper]=amppos;
                 phvel=phtmp[iper];
                 phtmp[iper]=phtmp[jper];
                 phtmp[jper]=phvel;
              }
           }
        ampst[nost]=(amptmp[0]/pertmp[0]+amptmp[1]/pertmp[1])/(1./pertmp[0]+1./pertmp[1]);
        phst[nost]=(phtmp[0]/pertmp[0]+phtmp[1]/pertmp[1])/(1./pertmp[0]+1./pertmp[1]);
        calc_dist( clatilst[i], clonglst[i], olati[nost], olong[nost], &dist );
        time[nost]=dist/phst[nost];
        if(dist<100) { csamp += ampst[nost]/pow(e_n,-1e-3*dist); nstmp++; }
//        calc_azimuth( latilst[i], longlst[i], olati[nost], olong[nost], &alpha[nost]);
        nost+=1;
/*        deg_l=(int)alpha[nost]; deg_h=deg_l+1;
        if(deg_h==360) deg_h=0;
        amp_n[deg_l][3]+=1; amp_n[deg_h][2]+=1;
        deg_l-=1; deg_h+=1;
        if(deg_l==-1) deg_l=359;
        if(deg_h==360) deg_h=0;
        amp_n[deg_l][3]+=1; amp_n[deg_h][2]+=1;

        d_alpha=alpha-deg_l+0.5;
        amp_n[deg_l][0]+=amppos/pow(e_n,-1e-3*dist)/d_alpha;
        amp_n[deg_l][1]+=1./d_alpha;
        amp_n[deg_l][3]+=1;
        d_alpha=deg_h-alpha+0.5;
        amp_n[deg_h][0]+=amppos/pow(e_n,-1e-3*dist)/d_alpha;
        amp_n[deg_h][1]+=1./d_alpha;
        amp_n[deg_h][2]+=1;  */
/*        if(deg_l==204) cout<<"ampsta: "<<amppos<<"  amp_h: "<<amp_n[deg_h][0]<<"  weight: "<<amp_n[deg_h][1]<<"  num: "<<amp_n[deg_h][2]<<endl;
        if(deg_l==205) cout<<"ampsta: "<<amppos<<"  amp_l: "<<amp_n[deg_l][0]<<"  weight: "<<amp_n[deg_l][1]<<"  num: "<<amp_n[deg_l][3]<<endl; */
     }
cout<<nost<<endl;
     if(nost==0) continue;
     if(nstmp==0) csamp = ampst[0]*3;
     else csamp /= nstmp;

/*
     for(k=0;k<20;k++)
        for(j=0;j<nost;j++) {
           if(alpha[j]*10<k+20 && alpha[j]*10>=k) flagh[k]=1;
           else if(alpha[j]*10<=k || alpha[j]*10>k+3580) flagl[k]=1;
        }
     for(k=20;k<3581;k++)
        for(j=0;j<nost;j++) {
           if(alpha[j]*10<k+20 && alpha[j]*10>=k) flagh[k]=1;
           else if(alpha[j]*10<=k && alpha[j]*10>k-20) flagl[k]=1;
        }
     for(k=3581;k<3600;k++)
        for(j=0;j<nost;j++) {
           if(alpha[j]*10<k-3580 || alpha[j]*10>=k) flagh[k]=1;
           else if(alpha[j]*10<=k && alpha[j]*10>k-20) flagl[k]=1;
        }
     for(k=0;k<3600;k++) flag[k]=flagh[k]*flagl[k];
*/

     sprintf(outname,"Ph_Amp_Map_%ssec/%s_center_ph_amp_map",argv[4], cstalst[i]);
     if((fout = fopen(outname, "w"))==NULL) {
       printf("Cannot open file %s to write\n",outname);
       continue;
     }
     fprintf(fout,"%9.5f %9.5f 000.00000 3.60000 %.5g %6s\n", clonglst[i], clatilst[i], csamp*1e8, cstalst[i] );
     for(j=0;j<nost;j++) {
//cout<<olong[j]<<"  "<<olati[j]<<"  "<<ampst[j]<<endl;
        fprintf(fout,"%9.5f %9.5f %10.5f %7.5f %.5g %6s\n", olong[j], olati[j], time[j], phst[j], ampst[j]*1e8, osta[nost] );
     }
     fclose(fout);

     sprintf(buff,"minmax %s | awk '{print $5,$6}' | sed s/'<'/''/g | sed s/'>'/''/g | sed s/'\\/'/' '/g | awk '{printf \"-R%%.0f/%%.0f/%%.0f/%%.0f\\n\",$1,$2,$3,$4}' > region_temp\0",outname);
     system(buff);
     sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel_positive %s region_temp 0.1\0", outname);
     //cout<<buff<<endl;
     system(buff);
     //exit(0);
     sprintf(buff,"mv %s.HD Ph_Amp_Map_%ssec/%s_center_ph_map.HD\0",outname,argv[4], cstalst[i]);
     system(buff);
     sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_ASN_am %s region_temp\0", outname);
     system(buff);
     sprintf(buff,"mv %s_am.HD Ph_Amp_Map_%ssec/%s_center_am_map.HD\0",outname,argv[4], cstalst[i]);
     system(buff);
     sprintf(buff,"mv %s_am Ph_Amp_Map_%ssec/%s_center_am_map\0",outname,argv[4], cstalst[i]);
     system(buff);
/*
     sprintf(filename,"%s.HD\0",outname);
     if((fin = fopen(filename, "r"))==NULL) {
       printf("Cannot open file %s to read\n",filename);
       continue;
     }
     for(j=0;;j++) {
        if( fgets(buff, 300, fin) == NULL ) break;
        if((sscanf(buff,"%f %f %f", &longtmp, &latitmp, &amppos))!=3) {
                 cout<<"Wrong format! Stopped: "<<filename<<endl;
                 break;
              }
        calc_dist( latilst[i], longlst[i], latitmp, longtmp, &dist );
        if(dist<100 || dist>500) continue;
        calc_azimuth( latilst[i], longlst[i], latitmp, longtmp, &alphatmp );
        alphatmp=(int)(alphatmp*10+0.5);
        if(alphatmp==3600) alphatmp=0;
        if(flag[(int)alphatmp]==0) continue;
        amp_n[(int)alphatmp][0]+=amppos;
        amp_n[(int)alphatmp][1]+=1;
     }
     fclose(fin);

     sprintf(outname,"Amp_Azimuth_%ssec/%s_amp_azimuth",argv[4],stalst[i]);
     if((fout = fopen(outname, "w"))==NULL) {
       printf("Cannot open file %s to write\n",outname);
       continue;
     }
     for(k=0;k<3600;k++) {
 //       cout<<"amp: "<<amp_n[k][0]<<"  num: "<<amp_n[k][1]<<endl;
        if(flag[k]==1 && amp_n[k][1]!=0) fprintf(fout,"%.1f   %.4g\n", k/10., amp_n[k][0]/amp_n[k][1]);
     }
     fclose(fout);
*/
     cout<<"The "<<i+1<<"th station "<<cstalst[i]<<" completed..."<<endl;
  }
  return 1;
}
