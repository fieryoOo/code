#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
using namespace std;


int STACK_dayflag_csta ( char *stalst, char *dirlst, char *outdir );

void hex2bin( char* hex, short* bin );

int calc_azimuth(double lati1, double long1, double lati2, double long2, double *alpha1);

int calc_dist(double lati1, double long1, double lati2, double long2, double *dist);

int Combiner (char *csta, char *dir, float *perlst, int nper, double *slong, double *slati, char **stalst, double **trvt, int *nsta);

int calc_grdt(char *file, double *slong, double *slati, double *sdat, double *grdt, double *azi, int *nsta);

//int calc_grdt_amp(char *stalist, char*csta, char *fflag, char *dirlist, double *slong, double *slati, double *sdat, double *grdt, double *azi, int *nsta);
int calc_lplc(char *csta, char *dir, float *perlst, int nper, double *slong, double *slati, double **dat_flag, double **grdt, double **azi, double **lplc, int *nsta);

int local_amp_grdt(char*csta, char *fsta, char *fflag, char *dirlist, float *perlst, int nper, double *slong, double *slati, char **stalst, double **sdat, double **grdt, double **azi, int *nsta);

#define NSTA 2000
int main (int argc, char *argv[])
{
   if(argc!=5){
      cout<<"usage:csta_map_gradient [csta_lst] [dir_lst] [All_Stack_dir] [per_lst]"<<endl;
      return 0;
   }
   FILE *ff,*ff2,*ff3;
   char fname[100], fname2[100], buff[300], *tsta[NSTA], **cstalst, *tmpc;
//   double grdt[NSTA], azit[NSTA], lplc[NSTA], amp[NSTA], grda[NSTA], azia[NSTA];
   double dec, prp;
   double tlong[NSTA], tlati[NSTA], along[NSTA], alati[NSTA];
   float *clonglst, *clatilst, per[50];
   int i, j, k, nper, ncsta, nsta, nstt;


   if((ff=fopen(argv[1],"r")) == NULL){
      cout<<"Can't open csta file: "<<argv[1]<<endl;
      return 0;
   }
   cstalst = NULL; clonglst = NULL; clatilst = NULL;
   for(ncsta=0;;ncsta++){
      cstalst=(char **) realloc (cstalst, (ncsta+1) * sizeof(char *));
      cstalst[ncsta] = (char *) malloc (10 * sizeof(char));
      clonglst = (float *) realloc (clonglst, (ncsta+1) * sizeof(float));
      clatilst = (float *) realloc (clatilst, (ncsta+1) * sizeof(float));
      if((fscanf(ff,"%s %f %f", cstalst[ncsta], &clonglst[ncsta], &clatilst[ncsta]))!=3) break;
      if(clonglst[ncsta]<0) clonglst[ncsta]+=360;
   }
   fclose(ff);
   for(i=0;i<NSTA;i++) tsta[i] = (char *) malloc (10 * sizeof(char));

   if((ff=fopen(argv[4],"r")) == NULL) {
      cout<<"Can't open per list file: "<<argv[4]<<endl;
      return 0;
   }
   for(nper=0;;nper++)
      if(fscanf(ff,"%f",&per[nper]) != 1)break;

   double *tvt[nper], *grdt[nper], *azit[nper], *lplc[nper];
   double *amp[nper], *grda[nper], *azia[nper];
   for(i=0;i<nper;i++) {
      tvt[i] = (double *) malloc (NSTA * sizeof(double));
      grdt[i] = (double *) malloc (NSTA * sizeof(double));
      azit[i] = (double *) malloc (NSTA * sizeof(double));
      lplc[i] = (double *) malloc (NSTA * sizeof(double));
      amp[i] = (double *) malloc (NSTA * sizeof(double));
      grda[i] = (double *) malloc (NSTA * sizeof(double));
      azia[i] = (double *) malloc (NSTA * sizeof(double));
   }


   for(i=0;i<ncsta;i++){
      cout<<"Working on the "<<i+1<<"th center station: "<<cstalst[i]<<endl;
//      sprintf(fname,"%s/Ph_Amp_Map_%ssec/%s_center_ph_amp_map_v2\0", argv[3], argv[4], argv[1]);
      cout<<"Extracting station and travel time information for all periods..."<<endl;
      for(j=0;j<nper;j++)for(k=0;k<NSTA;k++) tvt[j][k] = -1;
      Combiner (cstalst[i], argv[3], per, nper, tlong, tlati, tsta, tvt, &nstt);
      //double dec[nstt];
      cout<<"Computing travel time laplacian for all periods..."<<endl;
      calc_lplc(cstalst[i], argv[3], per, nper, tlong, tlati, tvt, grdt, azit, lplc, &nstt);
      //cout<<"nper: "<<nper<<" nsta: "<<nstt<<endl;
      //for(j=0;j<nstt;j++)cout<<tlong[j]<<" "<<tlati[j]<<" "<<tsta[j]<<" "<<tvt[0][j]<<endl;
      cout<<"Computing local amp gradient for all periods..."<<endl;
      sprintf(fname,"%s/Cor_dayflag.lst",argv[3]);
      sprintf(fname2,"%s/station.lst\0",argv[3]);
      local_amp_grdt(cstalst[i], fname2, fname, argv[2], per, nper, tlong, tlati, tsta, amp, grda, azia, &nstt);
//      for(j=0;j<nstt;j++) cout<<"ista: "<<j<<" "<<amp[0][j]<<" "<<grda[0][j]<<" "<<azia[0][j]<<endl;
      for(j=0;j<nper;j++) {
         sprintf(fname2,"%s_%.0fsec\0",cstalst[i],per[j]);
         sprintf(buff,"mkdir -p %s",fname2);
         system(buff);

         sprintf(fname,"%s/%s_prp_crd_%.0f",fname2,cstalst[i],per[j]);
         ff=fopen(fname,"w");
         sprintf(fname,"%s/%s_dec_%.0f",fname2,cstalst[i],per[j]);
         ff2=fopen(fname,"w");
         sprintf(fname,"%s/%s_prp_%.0f",fname2,cstalst[i],per[j]);
         ff3=fopen(fname,"w");
         for(k=0;k<nstt;k++) {
            if(azit[j][k]==-1 || azia[j][k]==-1) continue;
            dec=-grda[j][k]*cos(azia[j][k]-azit[j][k])/amp[j][k];
            prp=-0.5/grdt[j][k]*lplc[j][k];
            fprintf(ff,"%8.4f  %8.4f  %8g\n",tlong[k],tlati[k],dec+prp);
            fprintf(ff2,"%8.4f  %8.4f  %8g\n",tlong[k],tlati[k],dec);
            fprintf(ff3,"%8.4f  %8.4f  %8g\n",tlong[k],tlati[k],prp);
         }
         fclose(ff); fclose(ff2); fclose(ff3);

         sprintf(fname,"%s/%s_amp_%.0f",fname2,cstalst[i],per[j]);
         ff=fopen(fname,"w");
         sprintf(fname,"%s/%s_trvt_%.0f",fname2,cstalst[i],per[j]);
         ff2=fopen(fname,"w");
         sprintf(fname,"%s/%s_vel_%.0f",fname2,cstalst[i],per[j]);
         ff3=fopen(fname,"w");
         for(k=0;k<nstt;k++) {
            if(tvt[j][k]==-1) continue;
            fprintf(ff2,"%8.4f  %8.4f  %8g\n",tlong[k],tlati[k],tvt[j][k]);
         }
         for(k=0;k<nstt;k++) {
            if(azit[j][k]==-1) continue;
            fprintf(ff3,"%8.4f  %8.4f  %8g\n",tlong[k],tlati[k],1./grdt[j][k]);
         }
         for(k=0;k<nstt;k++) {
            if(amp[j][k]==-1) continue;
            fprintf(ff,"%8.4f  %8.4f  %8g\n",tlong[k],tlati[k],amp[j][k]);
         }
         fclose(ff); fclose(ff2); fclose(ff3);

      }
   }
return 1;
}
