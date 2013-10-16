#define MAIN
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
#include <string.h>
#include <iostream>
//#include <iomanip.h>
#include <math.h>
#include "/home/tianye/code/Programs/head/mysac64.h"
#define fMax 2000
using namespace std;


void fRemove (char *fname);

char * List(char *dir, const char *pattern, int type, int *nfile);

void FDivide (double f1, double f2, double f3, double f4, double dt, int n, float *seis_in, float *seis_out, double *freq, double *amp, double *pha, int nf);

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD) {
   FILE *fsac;
   if((fsac = fopen(fname, "r"))==NULL) return NULL;
   if ( !SHD ) SHD = &SAC_HEADER;
   fread(SHD,sizeof(SAC_HD),1,fsac);
   *sig = (float *) malloc (SHD->npts * sizeof(float));
   fread(*sig,sizeof(float),SHD->npts,fsac);
   fclose (fsac);

   /*-------------  calcule de t0  ----------------*/
   {
        int eh, em ,i;
        float fes;
        char koo[9];

        for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
        koo[8] = 0;

        SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
         SHD->nzsec + SHD->nzmsec*.001;

        sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

        SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}
   return SHD;
}

void write_sac (char *fname, float *sig, SAC_HD *SHD) {
   FILE *fsac;
   if( (fsac = fopen(fname, "wb"))==NULL ) {
      cerr<<"Cannot open file "<<fname<<endl;
      return;
   }
   if ( !SHD ) SHD = &SAC_HEADER;
   SHD->iftype = (int)ITIME;
   SHD->leven = (int)TRUE;
   SHD->lovrok = (int)TRUE;
   SHD->internal4 = 6L;

  /*+++++++++++++++++++++++++++++++++++++++++*/
   SHD->depmin = sig[0];
   SHD->depmax = sig[0];
   int i;
   for ( i = 0; i < SHD->npts ; i++ )
   {
    if ( SHD->depmin > sig[i] ) SHD->depmin = sig[i];
    else if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }

         fwrite(SHD,sizeof(SAC_HD),1,fsac);

         fwrite(sig,sizeof(float),SHD->npts,fsac);


        fclose (fsac);
}

void RTrend(float *sig, SAC_HD *shd) {
   // fit a*x+b
   int i, npts = shd->npts;
   float X = 0., Y = 0., X2 = 0., Y2 = 0., XY = 0.;
   for(i=0;i<npts;i++) {
      X += i;
      Y += sig[i];
      X2 += i*i;
      Y2 += sig[i]*sig[i];
      XY += i*sig[i];
   }
   float a = (npts*XY-X*Y)/(npts*X2-X*X);
   float b = (-X*XY+X2*Y)/(npts*X2-X*X);
   // correct sig and DEPMEN
   float mean = 0., max = sig[0], min = sig[0];
   float shift = b;
   for(i=0;i<npts;i++,shift+=a) {
      sig[i] -= shift;
      mean += sig[i];
      if ( min > sig[i] ) min = sig[i];
      else if ( max < sig[i] ) max = sig[i];
   }
   shd->depmin = min;
   shd->depmax = max;
   shd->depmen = mean / npts;
}

int TransferEvr(char *fname, char *fresp, char *erexe, float perl, float perh, float **sig, SAC_HD *sd) {
   // read in sac file
   SAC_HD shd;
   if( (read_sac(fname, sig, &shd))==NULL ) {
      fprintf(stderr, "*** Warning: cannot read sac file %s ***", fname);
      return 0;
   }
   *sd = shd;
   // running evalresp
   int nf = 100;
   char buff[300], sta[8], ch[8], net[8];
   float f2 = 1./perh, f1 = f2*0.9, f3 = 1./perl, f4 = f3*1.1;
   sscanf(sd->kstnm, "%s", sta);
   sscanf(sd->kcmpnm, "%s", ch);
   sscanf(sd->knetwk, "%s", net);
   sprintf(buff, "%s %s %s %4d %3d %f %f %d -f %s -v >& /dev/null", erexe, sta, ch, sd->nzyear, sd->nzjday, f1, f4, nf, fresp);
   system(buff);
   char nameam[50], nameph[50];
   sprintf(nameam, "AMP.%s.%s.*.%s", net, sta, ch);
   sprintf(nameph,"PHASE.%s.%s.*.%s", net, sta, ch);
   // find am file
   FILE *fam = NULL, *fph = NULL;
   int nlist;
   char *list = List(".", nameam, 0, &nlist);
   if( nlist!=1 ) {
      cerr<<"Error: "<<nlist<<" AMP file(s) found!"<<endl;
      exit(0);
   }
   sscanf(list, "%s", nameam);
   free(list);
   if( (fam = fopen(nameam, "r")) == NULL ) {
      cerr<<"Cannot open file "<<nameam<<endl;
      exit(0);
   }
   // find ph file
   list = List(".", nameph, 0, &nlist);
   if( nlist!=1 ) {
      cerr<<"Error: "<<nlist<<" PHASE file(s) found!"<<endl;
      exit(0);
   }
   sscanf(list, "%s", nameph);
   free(list);
   if( (fph = fopen(nameph, "r")) == NULL ) {
      cerr<<"Cannot open file "<<nameph<<endl;
      exit(0);
   }
   // read in am and ph data
   double pi=4*atan(1.0), pio180=pi/180.;
   double freq[nf], dtmp, amp[nf], pha[nf];
   int i = 0;
   while(i<nf) {
      if(fgets(buff, 300, fam)==NULL) break;
      sscanf(buff, "%lf %lf", &freq[i], &amp[i]);
      if(fgets(buff, 300, fph)==NULL) break;
      sscanf(buff, "%lf %lf", &dtmp, &pha[i]);
      if(dtmp!=freq[i]) {
	 cerr<<"incompatible AMP - PHASE pair!"<<endl;
	 continue;
      }
      amp[i] *= 0.000000001;
      pha[i] *= pio180;
      i++;
   }
   fclose(fam); fclose(fph);
   fRemove(nameam); fRemove(nameph);
   // remove trend ( and mean )
   RTrend(*sig, &shd);
   // run rmresponse
   FDivide (f1, f2, f3, f4, (double)shd.delta, shd.npts, *sig, *sig, freq, amp, pha, nf);
   return 1;
}

int main (int argc, char *argv[])
{
  if(argc != 6) {
    printf("Usage: %s [in_sac] [RESP_file] [out_sac] [perl] [perh]\n", argv[0]);
    exit(-1);
  }
  float *sig;
  SAC_HD sd;
  //double pi = 4*atan(1.0);

  float perl = atof(argv[4]), perh = atof(argv[5]);
  char erexe[100] = "/home/tianye/Software/evalresp/evalresp-3.3.3/evalresp";
  if( TransferEvr(argv[1], argv[2], erexe, perl, perh, &sig, &sd)!=1 ) {
     cerr<<"Transfer failed!"<<endl;
     exit(0);
  }

  write_sac (argv[3], sig, &sd );
  free(sig);

  return 1;  
}
