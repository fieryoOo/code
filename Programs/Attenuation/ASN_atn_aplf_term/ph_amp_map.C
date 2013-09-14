#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
using namespace std;

int calc_azimuth(double lati1, double long1, double lati2, double long2, double *alpha1);

int calc_dist(double lati1, double long1, double lati2, double long2, double *dist);


int ph_amp_map ( char *stafile, char *dayfile, float per, int daymin, int nrand)
{
  
  FILE *fin, *fin2, *fout;
  char filename[300], outname[300], buff[1000], buff2[300];
  char osta[5000][10], stalst[5000][10], csta[10], *tmp, *staname;
  char sta_pair[1000000][50], sta1[1000000][10], sta2[1000000][10];
  int i, j, k, iper, jper, nfile, nsta, ncsta, nost, nstmp;
  float tlength[1000000], ampst[5000], phst[5000], time[5000];
  float longlst[5000], latilst[5000], clong, clati;
  float per2, pertemp, perl, perh, pertmp[30], amptmp[30], phtmp[30];
  float period, amppos, snrpos, ampneg, snrneg, grvel, phvel, csamp; 
  float olong[5000], olati[5000], longtmp, latitmp, fstep=-0.1;
  double dist, alpha[5000], alphatmp, powe;
  double e_n=2.7182818284590451;
//  int flag[3600], flagl[3600]={0}, flagh[3600]={0};
//  if(argc != 6) 
//    {
//      printf("Usage: ASN_ph_amp_map [center_sta.lst] [sta.lst] [file_time-length.lst] [period] [time-length threshold]\n"); 
//      exit(-1);
//    }
  
//  per=atof(argv[4]);
  perl=1./exp(log(1./per)-fstep);
  perh=1./exp(log(1./per)+fstep);
  if((fin = fopen(dayfile, "r"))==NULL) {
    printf("Cannot open file list %s\n",dayfile);
    exit(1);
  }
  for(i=0;;i++) {
    if( (fgets(buff, 1000, fin)) == NULL ) break;
    sscanf(buff,"%s %f",&sta_pair[i], &tlength[i]);
//cout<<sta_pair[i]<<" "<<tlength[i]<<endl;
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

//  if((fin = fopen(argv[1], "r"))==NULL) {
//    printf("Cannot open center-station file %s\n",argv[1]);
//    exit(1);
//  }
//  for(i=0;;i++) {
//    if((fscanf(fin,"%s %f %f", cstalst[i], &clonglst[i], &clatilst[i]))!=3) break;
//    if(clonglst[i]<0) clonglst[i]+=360;
//  }
  
  ncsta = 1;
//  cout<<ncsta<<" stations read in. from center-station.lst"<<endl;
//  fclose(fin);

  if((fin = fopen(stafile, "r"))==NULL) {
    printf("Cannot open station file %s\n",stafile);
    exit(1);
  }
  fscanf(fin,"%s %f %f", csta, &clong, &clati);
  for(i=0;;i++) {
    if((fscanf(fin,"%s %f %f", stalst[i], &longlst[i], &latilst[i]))!=3) break;
    if(longlst[i]<0) longlst[i]+=360;
  }
  nsta = i;
  cout<<nsta<<" stations read in. from station.lst"<<endl;
  fclose(fin);

//  sprintf(buff,"mkdir -p Ph_Amp_Map_%.1fsec",per);
//  system(buff);
  for(i=0;i<ncsta;i++) {
//  for(i=1298;i<1299;i++) {
//     memset(amp_n,0,3600*2*sizeof(double));
     nost=0; nstmp=0; csamp=0;
     for(j=0;j<nfile;j++) {
        if(tlength[j]<daymin) continue;
        if(strcmp(csta,sta1[j])==0) {
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
           sprintf(filename,"STACK_%s/%s_amp_snr\0",stalst[0],sta_pair[j]);
           if((fin = fopen(filename, "r"))==NULL) {
//             printf("Amp file %s not found. Skipped\n",filename);
             continue;
           }
           sprintf(filename,"STACK_%s/%s_2_DISP.1\0",stalst[0],sta_pair[j]);
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
              if( snrpos>7 || (snrpos>5 && snrneg>7) ) {
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
        else if(strcmp(csta,sta2[j])==0) {
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
              if( snrneg>7 || (snrpos>7 && snrneg>5) ) {
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
        calc_dist( clati, clong, olati[nost], olong[nost], &dist );
        time[nost]=dist/phst[nost];
//        if(dist>atof(argv[5]) && dist<200) { cout<<ampst[nost]<<" "<<pow(e_n,-1e-3*dist)<<endl;if(ampst[nost]/pow(e_n,-1e-3*dist)>csamp) csamp = ampst[nost]/pow(e_n,-1e-3*dist); nstmp++; }
        if(dist>daymin && dist<200) { cout<<ampst[nost]<<" "<<pow(e_n,-1e-3*dist)<<endl; csamp += ampst[nost]/pow(e_n,-1e-3*dist); nstmp++; }
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
     if(nost==0) continue;
     if(nstmp==0) csamp = ampst[0]*2.5;
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
     sprintf(outname,"%d_ph_amp_map",nrand);
     if((fout = fopen(outname, "w"))==NULL) {
       printf("Cannot open file %s to write\n",outname);
       continue;
     }
     fprintf(fout,"%9.5f %9.5f  000.00000 3.50000 %6.5g %6s\n", clong, clati, csamp*1e8, csta );
     for(j=0;j<nost;j++) {
        fprintf(fout,"%9.5f %9.5f %10.5f %7.5f %6.5g %6s\n", olong[j], olati[j], time[j], phst[j], ampst[j]*1e8, osta[j] );
     }
     fclose(fout);

/*
     sprintf(buff,"minmax %s | awk '{print $5,$6}' | sed s/'<'/''/g | sed s/'>'/''/g | sed s/'\\/'/' '/g | awk '{printf \"-R%%.0f/%%.0f/%%.0f/%%.0f\\n\",$1,$2,$3,$4}' > region_temp\0",outname);
     system(buff);
     sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel_positive %s region_temp 0.1\0", outname);
     //cout<<buff<<endl;
     system(buff);
     //exit(0);
     sprintf(buff,"mv %s.HD %s_center_ph_map.HD\0",outname, csta);
     system(buff);
     sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_ASN_am %s region_temp\0", outname);
     system(buff);
*/
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
//     cout<<"The "<<i+1<<"th station "<<csta<<" completed..."<<endl;
  }
  return 1;
}
