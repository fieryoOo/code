#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;

#define BLKSIZE 2000

struct DATA {
   float val;
   float sigma;
};

struct DATA *data;

int readdata( char* fname ) {
   FILE *fin;
   if( (fin=fopen(fname,"r")) == NULL ) {
      cerr<<"Cannot access file "<<fname<<endl;
      exit(-1);
   }
   int nblk = 0, ndat;
   char buff[300];
   data = NULL;
   for(ndat=0;fgets(buff, 300, fin);) {
      if( nblk*BLKSIZE <= ndat ) data = (struct DATA *) realloc (data, (++nblk)*BLKSIZE * sizeof(struct DATA));
      if( sscanf(buff, "%f %f", &(data[ndat].val), &(data[ndat].sigma)) != 2 ) {
         cerr<<"Warning(main): format error: "<<buff<<endl;
         continue;
      }
      ndat++;
   }
   fclose(fin);
   return ndat;
}

float AvgSTD(struct DATA *data, int ndat, float *stdo) {
   int i;
   // first iteration
   float avg=0., std=0., V1=0., V2=0.;
   float ftmp, weight[ndat];
   for(i=0;i<ndat;i++) {
      weight[i] = 1. / sqrt(data[i].sigma);
      V1 += weight[i]; V2 += weight[i] * weight[i];
   }
   for(i=0; i<ndat; i++) avg += data[i].val * weight[i];
   avg /= V1;
   for(i=0;i<ndat;i++) {
      ftmp = data[i].val - avg;
      std += ftmp*ftmp*weight[i];
   }
   std = sqrt(std * V1 / (V1*V1-V2) );
   //cerr<<avg<<" "<<std<<endl;
   // second iteration
   V1=0., V2=0.;
   for(i=0;i<ndat;i++) {
      if( fabs(data[i].sigma - avg) > std*2 ) weight[i] = 0.;
      else { V1 += weight[i]; V2 += weight[i] * weight[i]; }
   }
   avg=0., std=0.;
   for(i=0; i<ndat; i++) avg += data[i].val * weight[i];
   avg /= V1;
   float nwstd=0.;
   for(i=0;i<ndat;i++) {
      ftmp = data[i].val - avg;
      std += ftmp*ftmp*weight[i];
      nwstd += ftmp*ftmp;
   }
   std = sqrt(std * V1 / (V1*V1-V2) );
   nwstd = sqrt(nwstd / (ndat-1));
   *stdo = std;
   //cerr<<avg<<" "<<std<<endl;
   return nwstd;
}

void KernelDensity(struct DATA *data, int ndat, float std, float h, char *outname) {
   float step = std/50.;
   int n5sig = int(5*std/step+0.5), nden = 2 * n5sig + 1;
   float amp, alpha, density[nden];
   float sigma, val, ftmp;
   float osqrt2pi = 1./sqrt(2.*4.*atan(1.));
   int i, iden, ib, ie;
   memset(density, 0, nden*sizeof(float));
   for(i=0;i<ndat;i++) {
      sigma = data[i].sigma * h;
      val = data[i].val;
      amp = osqrt2pi/sigma;
      //cerr<<i<<": "<<amp<<" "<<sigma<<" "<<density[n5sig]<<endl;
      alpha = 0.5/(sigma*sigma);
      ib = (int)floor((val-3.5*sigma)/step) + n5sig;
      ie = (int)ceil((val+3.5*sigma)/step) + n5sig + 1;
      if( ib < 0 ) ib = 0;
      if( ie > nden ) ie = nden;
      for(iden=ib; iden<ie; iden++) {
	 ftmp = (val - (iden-n5sig)*step);
	 density[iden] += amp*exp(-alpha*ftmp*ftmp);
      }
   }
   FILE *fout = fopen(outname, "w");
   for(iden=0;iden<nden;iden++) fprintf(fout, "%f %f\n", (iden-n5sig)*step, density[iden]/ndat);
   fclose(fout);
}

int main(int argc, char *argv[]) {
   if(argc != 2) {
      cerr<<"Usage: "<<argv[0]<<" [input_file (col1=value col2=uncertainty)]"<<endl;
      exit(-1);
   }
   // read in data
   int ndat = readdata(argv[1]);
   // compute bandwidth using Gaussian approximation
   float std, nwstd = AvgSTD(data, ndat, &std);
   float h = std*pow(1.3333333/ndat, 0.2);
   cout<<std<<" "<<nwstd<<" "<<h<<endl;
   // estimate kernel density and write to file
   char outname[100];
   sprintf(outname, "%s_kd", argv[1]);
   KernelDensity(data, ndat, std, h, outname);

   return 0;
}
