#define MAIN
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include "/home/tianye/code/Programs/head/mysac64.h"
#define fMax 2000
using namespace std;

void rmresponse_(int *n,double *dt,float *sei, double *freq, double *phase_res, double *amp_res, int *nf);

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
   fsac = fopen(fname, "wb");
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
    if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }

         fwrite(SHD,sizeof(SAC_HD),1,fsac);

         fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);


        fclose (fsac);
}


int main (int argc, char *argv[])
{
   if(argc != 4) {
      printf("Usage: do_cor [in_sac_list (file1 file2)] [time_shift] [out_file]\n");
      exit(-1);
   }

   float *sigc, *sigo;
   SAC_HD sdc, sdo;
   FILE *flst, *ff;
   char buff[300], sacc[150], saco[150];
   int i, ii, bc, bo, dnum, ncut = 0;
   int lag = 1000, nhalf = 0;
   float cor[2*lag*50+1], ftemp;
   
   if( (flst=fopen(argv[1], "r"))==NULL ) {
      cout<<"Cannot access file "<<argv[1]<<endl;
      exit(0);
   }
   for(i=0;i<2*lag*50+1;i++) cor[i] = 0.;
   if( fgets(buff, 300, flst)!=NULL ) {
     if( sscanf(buff, "%s %s", sacc, saco)!=2 ) {
        cout<<"File list format error!"<<endl;
        exit(0);
     }
     if( (read_sac(sacc, &sigc, &sdc))==NULL ) {
         fprintf(stderr, "*** Warning: cannot read sac file %s ***\n", sacc);
         exit(0);
     }
     if( (read_sac(saco, &sigo, &sdo))==NULL ) {
         fprintf(stderr, "*** Warning: cannot read sac file %s ***\n", saco);
         exit(0);
     }
     if(sdc.delta!=sdo.delta || sdc.npts!=sdc.npts) {
        cout<<"*** Warning: header file mismach ***"<<endl;
        exit(0);
     }
     if(nhalf == 0) nhalf = (int)floor(lag/sdc.delta);
     int nshift = (int)floor(atof(argv[2])/sdc.delta+0.5);
     //else if(nhalf != (int)floor(lag/sdc.delta)) continue;
     ff = fopen(argv[3], "w");
     if(nshift>=0) for(ii=0,bc=ncut,bo=nshift+ncut; bo<sdo.npts-ncut; bo++,bc++) {
        ftemp = (sigc[bc] / sigo[bo]);
	fprintf(ff, "%f %f %f %f\n", bc*0.2, sigc[bc], sigo[bo], ftemp);
        //ii++;
     }
     else for(ii=0,bc=ncut-nshift,bo=ncut; bc<sdc.npts-ncut; bo++,bc++) {
        ftemp = (sigc[bc] / sigo[bo]);
	fprintf(ff, "%f %f %f %f\n", bc*0.2, sigc[bc], sigo[bo], ftemp);
        //ii++;
     }
     fclose(ff);
   }
   fclose(flst);

   return 1;  
}
