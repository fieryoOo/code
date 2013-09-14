#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
using namespace std;

#define DMIN 60

void hex2bin( char* hex, short* bin );

int STACK_dayflag_csta ( char **staname, float *stalon, float *stalat, int nsta, char **dirlist, int ndir, char *outdir );

int ph_amp_map ( char *stafile, char *dayfile, float per, int daymin, int nrand);

int local_amp_grdt(char*csta, char *fsta, char *fflag, char *dirlist, float *perlst, int nper, double *slong, double *slati, char **stalst, double **sdat, double **grdt, double **azi, int *nsta)
{
   FILE *ff;
   char buff[300], buff2[300], stabuff[300], bufflg[1000], dirlst[100][100], *staname[100], *dirname[50];
   float clong, clati;
   char stnm1[10], stnm2[10], stal[2000][10];
   float longmin=360, longmax=0, latimin=90, latimax=-90, longtmp, latitmp;
   float stalon[100], stalat[100];
   double stalo[2000], stala[2000];
   double radius=6371.1391285, pi=4.0*atan(1.0), dgr=0.1, dattmp;
   double grdtx, grdty, gradtx, gradty, dis, alpha, weight;
   int i, j, ii, jj, kk, iii, jjj, ista, nstl, nrand, ndir, nSsta, nSdir, ntmpp;

//   sprintf(stalst[0],"%s",csta);
//   if((ff=fopen(ph_stalist,"r"))==NULL){
//     printf("Can't open %s to read!\n",ph_stalist);
//     exit(0);
//   }
   i=0;
   if(strcmp(stalst[0],csta) != 0){
//      for(i=1;i<*nsta;i++) if(strcmp(stalst[i],csta) == 0) {
//         clong=slong[i]; clati=slati[i];
//         slong[i]=slong[0]; slati[i]=slati[0];
//         sprintf(stalst[i],"%s\0",stalst[0]);
//         break;
//      }
//      if(i==*nsta) {
//         cout<<"Cannot find center station "<<csta<<" in stalst"<<endl;
      cout<<"The center station has to be the 1st one in the input station list. Stopped!"<<endl;
      exit(0);
//      }
   }
   else { clong=slong[0]; clati=slati[0]; }

   if((ff=fopen(fsta,"r")) == NULL) {
      cout<<"Can't open station list: "<<fsta<<endl;
      exit(0);
   }
   for(nstl=0;;nstl++){
      if(fgets(buff, 300, ff) == NULL) break;
      sscanf(buff,"%s %lf %lf", &stal[nstl], &stalo[nstl], &stala[nstl]);
      if(stalo[nstl] < 0) stalo[nstl] += 360;
   }
   fclose(ff);

/*
   flag=0;
   for(*nsta=1;;){
      if(fgets(buff, 300, ff) == NULL ) break;
      sscanf(buff,"%lf %lf %lf %lf %lf %s", &slong[*nsta], &slati[*nsta], &dattmp, &dattmp, &dattmp, &stalst[*nsta]);
      if(slong[*nsta] < 0) slong[*nsta] += 360;
      if(strcmp(csta,stalst[*nsta]) == 0) {
         clong=slong[*nsta]; clati=slati[*nsta];
         flag=1;
         continue;
      }
      *nsta=*nsta+1;
   }
   fclose(ff);
   if(flag==0) {
      cout<<"Cannot find center station "<<csta<<" in "<<ph_stalist<<endl;
      return 0;
   }
*/

   if((ff=fopen(dirlist,"r"))==NULL){
      printf("Can't open %s to read!\n",dirlist);
      exit(0);
   }
   for(ndir=0;;ndir++) {
      if(fgets(buff, 300, ff) == NULL ) break;
      sscanf(buff,"%s",&dirlst[ndir]);
   }

   char d_flag[nstl][900], day_flag[900];
   for(ista=1;ista<nstl;ista++) sprintf(d_flag[ista],"-1\0");
   if((ff=fopen(fflag,"r")) == NULL){
     printf("Can't open %s to read!\n",fflag);
     exit(0);
   }
   for(;;){
      if(fgets(bufflg, 1000, ff) == NULL ) break;
      strtok(bufflg,"\n");
      if(sscanf(bufflg,"%s %d %[^/n]", &buff, &ntmpp, &day_flag) != 3) {
         cout<<"Wrong format! Stopped: "<<bufflg<<endl;
         break;
      }
      sscanf(buff, "%[^_]_%[^_]_%[^.]", &buff2, &stnm1, &stnm2);
      if(strcmp(csta,stnm1) == 0) {
        for(i=0;i<nstl;i++) if(strcmp(stnm2,stal[i]) == 0) break;
        if(i == nstl) continue;
        sprintf(d_flag[i], "%s", day_flag);
      }
      else if(strcmp(csta,stnm2) == 0) {
        for(i=0;i<nstl;i++) if(strcmp(stnm1,stal[i]) == 0) break;
        if(i == nstl) continue;
        sprintf(d_flag[i], "%s", day_flag);
      }
   }
   fclose(ff);

/*
   int npts_long=int((longmax-longmin)/dgr+1);
   int npts_lati=int((latimax-latimin)/dgr+1);
   double dislong[npts_lati], dislati, distmp;
   double gradtr[npts_long][npts_lati], grdtmp, azig[npts_long][npts_lati], dat[npts_long][npts_lati];

   for(i=0;i<npts_long;i++) for(j=0;j<npts_lati;j++){
      gradtr[i][j] = 1e10;
      dat[i][j] = -1;
   }
   for(i=1;i<npts_lati-1;i++){
      dislati=atan(0.993277*tan((latimin+i*dgr)/180*pi))*180/pi;
      dislong[i]=radius*sin((90-dislati)/180*pi)*dgr/180*pi;
   }
   dislati=radius*dgr/180*pi;
*/

   int npts_long, npts_lati, imax, nend;
   double gradtr, grdtmp, dat[4][4];
   double dislong, dislati=radius*dgr/180*pi, distmp;
   char d_hex[9], *hex_tmp;
   short nmonth, nnmonth, nmax, ndata, nii, njj, ns11, ns12, ns21, ns22;
   short cd_bin[3000]={0}, d_bin[32]={0};
   short sta_marker[100], month_marker[ndir], nosta, month_flag[100][ndir], n_month[100];
   int cday[ndir], sday, osta;
   srand ( time(NULL) );
   nrand=rand();

   for(ista=0;ista<*nsta;ista++) for(j=0;j<nper;j++) sdat[j][ista] = -1.;

   for(ista=1;ista<*nsta;ista++){
      for(ii=0;ii<nper;ii++)azi[ii][ista]=-1;
      for(osta=0;osta<nstl;osta++) 
         if((strcmp(stalst[ista],stal[osta])) == 0) break;
      if(osta == nstl) continue;
      if(strcmp(d_flag[osta],"-1") == 0) continue;
      dislong=atan(0.993277*tan(slati[ista]/180*pi))*180/pi;
      dislong=radius*sin((90-dislong)/180*pi)*dgr/180*pi;
      latimin=slati[ista]-1.8; latimax=slati[ista]+1.8;
      longmax=200./(dislong/dgr);
      longmin=slong[ista]-longmax; longmax=slong[ista]+longmax;
/*
      sprintf(buff,"region_%d",nrand);
      ff=fopen(buff,"w");
      fprintf(ff, "-R%d/%d/%d/%d\0", (int)floor(longmin), (int)ceil(longmax), (int)floor(latimin), (int)ceil(latimax));
      fclose(ff);
*/
      sprintf(day_flag,"%s",d_flag[osta]);
      hex_tmp=strtok(day_flag," ");
      sprintf(d_hex,"%s",hex_tmp);
      hex2bin(d_hex,&cd_bin[0]);
      cday[0]=0;
      for(ii=1;ii<32;ii++) cday[0] += cd_bin[ii];
      if(cday[0]>15) month_flag[0][0]=1;
      else month_flag[0][0]=0;
      for(i=1;;i++){
         if((hex_tmp=strtok(NULL," ")) == NULL) break;
         sprintf(d_hex,"%s",hex_tmp);
         hex2bin(d_hex,&cd_bin[i*32]);
         cday[i]=0;
         for(ii=1;ii<32;ii++) cday[i] += cd_bin[i*32+ii];
         if(cday[i]>15) month_flag[0][i]=1;
         else month_flag[0][i]=0;
      }
      nmonth=i;
      if(nmonth != ndir) {
         cout<<"Number of months mismatch between dir.lst and dayflag.lst!"<<endl;
         exit(0);
      }
//      sprintf(buff,"%s_%d", stalst[ista], nrand);
//      ff=fopen(buff,"w");
      nosta=1;
      sta_marker[0] = osta;
      for(i=0;i<nstl;i++){
         //if(i==ista) continue;
         if(i == osta) continue;
         //if(slati[i]<latimin || slati[i]>latimax || slong[i]<longmin || slong[i]>longmax) continue;
         if(stala[i]<latimin || stala[i]>latimax || stalo[i]<longmin || stalo[i]>longmax) continue;
         sprintf(day_flag,"%s",d_flag[i]);
         sprintf(d_hex,"%s",strtok(day_flag," "));
         hex2bin(d_hex,d_bin);
         sday=0;
         for(j=1;j<32;j++) if(d_bin[j]==1) {
            if(cd_bin[j] == 1) sday += 1;
            else sday -= 1;
         }
         if((float)sday/(float)cday[0]>0.8) month_flag[nosta][0]=1;
         else month_flag[nosta][0]=0;
         for(ii=1;ii<nmonth;ii++){
            if((hex_tmp=strtok(NULL," ")) == NULL) break;
            if(month_flag[0][ii]==0) continue;
            sprintf(d_hex,"%s",hex_tmp);
            hex2bin(d_hex,d_bin);
            sday=0;
            for(j=1;j<32;j++) if(d_bin[j]==1) {
               if(cd_bin[ii*32+j] == 1) sday += 1;
               else sday -= 1;
            }
            if((float)sday/(float)cday[ii]>0.8) month_flag[nosta][ii]=1;
            else month_flag[nosta][ii]=0;
         }
         j=0;
         for(ii=0;ii<nmonth;ii++)
            if(month_flag[0][ii]==1 && month_flag[nosta][ii]==1) j += 1;
         if(j<3) continue;
         sta_marker[nosta] = i;
//         if(slati[i]<slati[ista] && slong[i]<slong[ista]) ns11++;
//         else if(slati[i]<slati[ista] && slong[i]>slong[ista]) ns12++;
//         else if(slati[i]>slati[ista] && slong[i]<slong[ista]) ns21++;
//         else ns22++;
//         sday=0;
//         fprintf(ff,"%lf %lf %lf\n", slong[i], slati[i], sdat[i]);
         nosta++;
      }
//      fclose(ff);
/*
      for(ii=0;ii<nosta;ii++) {
         n_month[ii]=0;
         for(jj=0;jj<nmonth;jj++) n_month[ii] += month_flag[ii][jj];
      }
      imax=1;
      for(jj=1;jj<nosta;jj++) if(n_month[imax]<n_month[jj]) imax=jj;
      if(n_month[0]<n_month[imax]){
         j=n_month[imax];
         n_month[imax]=n_month[0];
         n_month[0]=j;
         j=sta_marker[imax];
         sta_marker[imax]=sta_marker[0];
         sta_marker[0]=j;
         for(jj=0;jj<nmonth;jj++){
            j=month_flag[0][jj];
            month_flag[0][jj]=month_flag[imax][jj];
            month_flag[imax][jj]=j;
         }
      }
*/
/*/out:
      cout<<"     ";
      for(jj=0;jj<nmonth;jj++){
          if(month_flag[0][jj]==0)continue;
          printf("%2d",jj);
      }
      cout<<endl;
      for(ii=0;ii<nosta;ii++){
          printf("%-5s",stalst[sta_marker[ii]]);
          for(jj=0;jj<nmonth;jj++){
             if(month_flag[0][jj]==0)continue;
             printf("%2d",month_flag[ii][jj]);
          }
          cout<<endl;
      }
//:out*/

      for(ii=1;ii<nosta-1;ii++){
         imax=ii;
         for(kk=ii;kk<nosta;kk++){
            n_month[kk]=0;
            for(jj=0;jj<nmonth;jj++) {
               if(month_flag[0][jj]==0)continue;
               if(month_flag[kk][jj] == 1) {
                  if(month_flag[ii-1][jj]==1)n_month[kk] += 1;
                  else n_month[kk] -= 1;
               }
            }
            if(n_month[imax]<n_month[kk]) imax=kk;
         }
         if(n_month[ii]<n_month[imax]){
            //j=n_month[imax];
            //n_month[imax]=n_month[ii];
            //n_month[ii]=j;
            j=sta_marker[imax];
            sta_marker[imax]=sta_marker[ii];
            sta_marker[ii]=j;
            for(jj=0;jj<nmonth;jj++){
               j=month_flag[ii][jj];
               month_flag[ii][jj]=month_flag[imax][jj];
               month_flag[imax][jj]=j;
            }
         }
      }

      for(ii=0;ii<nmonth;ii++) {
         if(month_flag[0][ii]==0)continue;
         for(jj=2;jj<nosta;jj++) if(month_flag[jj][ii]==1)
            month_flag[jj][ii]+=month_flag[jj-1][ii];
      }
/*/out:
      cout<<"     ";
      for(jj=0;jj<nmonth;jj++){
          if(month_flag[0][jj]==0)continue;
          printf("%2d",jj);
      }
      cout<<endl;
      for(ii=0;ii<nosta;ii++){
          printf("%-5s",stalst[sta_marker[ii]]);
          for(jj=0;jj<nmonth;jj++){
             if(month_flag[0][jj]==0)continue;
             printf("%2d",month_flag[ii][jj]);
          }
          cout<<endl;
      }
//:out*/


      for(ii=1;ii<nosta;ii++) {
         for(kk=6;kk<ii;kk++) {
            ns11=0; ns12=0; ns21=0; ns22=0;
            for(jj=ii;jj>=ii-kk;jj--){
               if(stala[sta_marker[jj]]<stala[osta] && stalo[sta_marker[jj]]<stalo[osta]) ns11++;
               else if(stala[sta_marker[jj]]<stala[osta] && stalo[sta_marker[jj]]>stalo[osta]) ns12++;
               else if(stala[sta_marker[jj]]>stala[osta] && stalo[sta_marker[jj]]<stalo[osta]) ns21++;
               else ns22++;
            }
            if(ns11>0 && ns12>0 && ns21>0 && ns22>0) break;
         }
         for(jj=0;jj<nmonth;jj++)
            if(month_flag[ii][jj]<kk+1) month_flag[ii][jj]=0;
      }

      nmax=0;
      for(ii=7;ii<nosta;ii++) {
         for(jj=0;jj<nmonth;jj++) {
            if(month_flag[0][jj]==0)continue;
            if(month_flag[ii][jj]==0) continue;
            nnmonth=0;
            for(kk=0;kk<nmonth;kk++)
               if(month_flag[0][kk]!=0 && month_flag[ii][kk]>=month_flag[ii][jj]) nnmonth++;
            if(nnmonth<3) continue;
            ndata=month_flag[ii][jj]*nnmonth;
            if(nmax<ndata) { nmax=ndata; nii=ii; njj=jj; }
         }
      }
//cout<<"month*sta: "<<nmax<<" sta#: "<<nii<<" month#: "<<njj<<endl;
/*out:
      cout<<"     ";
      for(jj=0;jj<nmonth;jj++){
          if(month_flag[0][jj]==0)continue;
          printf("%2d",jj);
      }
      cout<<endl;
      for(ii=0;ii<nosta;ii++){
          printf("%-5s",stalst[sta_marker[ii]]);
          for(jj=0;jj<nmonth;jj++){
             if(month_flag[0][jj]==0)continue;
             printf("%2d",month_flag[ii][jj]);
          }
          cout<<endl;
      }
//:out*/

      if(nmax==0) continue;
      staname[0]=&csta[0]; stalon[0]=clong; stalat[0]=clati;
      staname[1]=&stalst[ista][0]; stalon[1]=slong[ista]; stalat[1]=slati[ista];

      nSsta=2;
      for(ii=nii-month_flag[nii][njj]+1;ii<=nii;ii++){
         staname[nSsta] = &stal[sta_marker[ii]][0];
         stalon[nSsta] = stalo[sta_marker[ii]];
         stalat[nSsta] = stala[sta_marker[ii]];
         nSsta++;
      }
      nSdir=0;
      for(jj=0;jj<nmonth;jj++){
         if(month_flag[0][jj]==0)continue;
         if(month_flag[nii][jj]<month_flag[nii][njj])continue;
         dirname[nSdir] = &dirlst[jj][0];
         nSdir++;
      }


      sprintf(buff,"STACK_%s",stalst[ista]);
      STACK_dayflag_csta ( staname, stalon, stalat, nSsta, dirname, nSdir, buff );
      sprintf(stabuff,"STACK_%s/station.lst",stalst[ista]);
      ff=fopen(stabuff,"w");
      for(jj=0;jj<nSsta;jj++) fprintf(ff,"%s\t%lf\t%lf\n",staname[jj],stalon[jj],stalat[jj]);
      fclose(ff);
      longmin=floor(longmin); longmax=ceil(longmax);
      latimin=floor(latimin); latimax=ceil(latimax);
      sprintf(buff,"region_%d",nrand);
      ff=fopen(buff,"w");
      fprintf(ff, "-R%d/%d/%d/%d\n", (int)longmin, (int)longmax, (int)latimin, (int)latimax);
      fclose(ff);

      sprintf(buff,"/home/tianye/code/Programs/Attenuation/ASN_atn_aplf_term/FTAN.csh STACK_%s/%s\0", stalst[ista], csta);
      system(buff);

      for(jj=0;jj<nper;jj++) {
cout<<"di "<<jj<<" ge period"<<endl;
         sprintf(buff, "STACK_%s/Cor_dayflag.lst\0", stalst[ista]);
         ph_amp_map(stabuff, buff, perlst[jj], DMIN, nrand);
         sprintf(buff,"/home/tianye/code/Programs/Attenuation/ASN_atn_aplf_term/correct_2pi_curvature %d %f\0", nrand, perlst[jj]);
         system(buff);
//         sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_travel_positive %d_ph_amp_map region_%d 0.1\0", nrand, nrand);
//         system(buff);
//         sprintf(buff,"mv %d_ph_amp_map.HD %d_ph_map.HD\0",nrand, nrand);
//         system(buff);
//         sprintf(buff,"/home/tianye/code/Script/GMT/C_plot_ASN_am %d_ph_amp_map region_%d\0", nrand,nrand);
//         system(buff);

         sprintf(buff2,"%d_am_map_v2",nrand);
         if((ff=fopen(buff2,"r"))==NULL){
            continue;
         }
         for(;;){
            if( fgets(buff2, 300, ff) == NULL ) break;
            sscanf(buff2,"%f %f %lf %s", &longtmp, &latitmp, &dattmp, buff);
            if(strcmp(buff, stalst[ista]) == 0) { 
               sdat[jj][ista]=dattmp;
               break; 
            }
         }
         fclose(ff);
         if(sdat[jj][ista] == -1) continue;
         sprintf(buff2,"%d_am_map_v2.HD",nrand);
         if((ff=fopen(buff2,"r"))==NULL){
            continue;
         }
         for(;;){
            if( fgets(buff2, 300, ff) == NULL ) break;
            sscanf(buff2,"%f %f %lf", &longtmp, &latitmp, &dattmp);
            i=(int)floor((longtmp-slong[ista])/dgr)+2;
            j=(int)floor((latitmp-slati[ista])/dgr)+2;
            if(i<0||i>3||j<0||j>3) continue;
            //if(i<0||i>=npts_long||j<0||j>=npts_lati) continue;
            dat[i][j]=dattmp;
         }
         fclose(ff);
         gradtx = 0; gradty = 0; ntmpp=0; weight=0;
         for(i=1;i<3;i++)
            for(j=1;j<3;j++){
               distmp=sqrt(dislong*dislong+dislati*dislati);
               grdtx=(dat[i+1][j-1]-dat[i-1][j+1])/2.0/distmp;
               grdty=(dat[i+1][j+1]-dat[i-1][j-1])/2.0/distmp;
               alpha=2*atan(dislati/dislong);
               alpha=atan((grdtx/grdty-cos(alpha))/sin(alpha));
               grdtmp=fabs(grdty/cos(alpha));
               grdtx=(dat[i+1][j]-dat[i-1][j])/2.0/dislong;
               grdty=(dat[i][j+1]-dat[i][j-1])/2.0/dislati;
               gradtr=sqrt(grdtx*grdtx+grdty*grdty);
               if(fabs((gradtr-grdtmp)/gradtr)>0.15) continue;
               dis=sqrt(pow((i-1-slong[ista]+(int)slong[ista])*dislong,2)+pow((j-1-slati[ista]+(int)slati[ista])*dislati,2))+1e-20;
               gradtx += grdtx/dis; gradty += grdty/dis;
               //grdt[ista]+=gradtr[ntmpp]/dis;
               //azi[ista]+=azig[ntmpp]/dis;
               weight+=1/dis;
               ntmpp++;
            }
         if(ntmpp<2) { azi[jj][ista]=-1; continue; }
         grdt[jj][ista]=sqrt(gradtx*gradtx+gradty*gradty)/weight;
         azi[jj][ista]=pi/2-atan2(gradty,gradtx);
         if(azi[jj][ista]<0-1e-8)azi[jj][ista]+=2*pi;
     }
     sprintf(buff,"rm -f region_%d",nrand);
     system(buff);
     sprintf(buff,"rm -f %d_*",nrand);
     system(buff);
     sprintf(buff,"rm -rf STACK_%s\0",stalst[ista]);
     system(buff);
  }

  return 1;
}

