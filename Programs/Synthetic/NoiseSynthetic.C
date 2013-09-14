//The fk synthetic package called from this program doesn't seem to be thread-safe. Debug needed!

#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include <errno.h>
#include "mysac64.h"
#include "DisAzi.h"
using namespace std;

#define NBLKsrc 500
#define NBLKsta 10
#define PIO180 0.017453292519943295

#define NTHRDS 4

char buff[300];
int nsrc, nsta, npts = 86401, dep = 5;
int segnpt = 500, segnum = (int)ceil((float)npts/segnpt);
float dt = 1.;
struct SRC *src = NULL;
struct STA *sta = NULL;
struct PAIR **pairs = NULL;
pthread_mutex_t **lock;

struct SRC {
   float lon, lat;
   float amp, foccur;
};

struct STA {
   float lon, lat;
   char name[50];
   float *sig;
   SAC_HD shd;
};

struct PAIR {
   float dis, azi;
   float *sigx, *sigy, *sigz;
   SAC_HD shdx, shdy, shdz;
};

//struct iRNG {
//   int b, e;
//};


SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD) {
   FILE *fsac;
   if((fsac = fopen(fname, "rb"))==NULL) return NULL;
//   if ( !SHD ) SHD = &SAC_HEADER;
   fread(SHD,sizeof(SAC_HD),1,fsac);
   //if( *sig == NULL) 
   *sig = (float *) malloc (SHD->npts * sizeof(float));
   //else *sig = (float *) realloc (*sig, SHD->npts * sizeof(float));
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
   if( (fsac = fopen(fname, "wb"))==NULL ) cout<<"Cannot open file "<<fname<<endl;
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

uint64_t ClockGetTime()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
}
int RandSRC(float prob, float *stk, float *dip) {
   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
//   boost::normal_distribution<> normal(0,1);
//   boost::variate_generator< RNGType, boost::normal_distribution<> > Gauss(gener, normal);
   boost::uniform_real<> uniform(0.,1.);
   boost::variate_generator< RNGType, boost::uniform_real<> > Rand(gener, uniform);
   if(Rand()>prob) return 0;
   *stk = 360.*Rand();
   *dip = 90.*(2.*Rand()-1.);
   return 1;
}

int ClipAddSac(int isrc, int ista, float stk, float dip, float tshift, float **sig2, SAC_HD *sd) {
   struct PAIR pair = pairs[isrc][ista];
   float amplf = src[isrc].amp;

   stk *= PIO180; dip *= PIO180;
   int i, nb, ne;
   float ux = cos(dip)*sin(stk)*amplf, uy = cos(dip)*cos(stk)*amplf, uz = sin(dip)*amplf;
   long double fra1, fra2;
   nb = (int)ceil((pair.shdz.b+tshift)/pair.shdz.delta);
   fra1 = nb-(pair.shdz.b+tshift)/pair.shdz.delta;
   fra2 = 1.-fra1;
//   if( npts<nb+shd.npts-1 ){
//      cout<<"Increase npts!"<<endl;
//      exit(0);
//   }
   //compute syn signal and store in local array
   float sigtmp[pair.shdz.npts];
   i = 0;
   if(nb<0) i = -nb;
   if(fra1==0)
      for(;i<pair.shdz.npts-1;i++) sigtmp[i] = ux*pair.sigx[i]+uy*pair.sigy[i]+uz*pair.sigz[i];
   else
      for(;i<pair.shdz.npts-1;i++) sigtmp[i] = ux*(pair.sigx[i]*fra2 + pair.sigx[i+1]*fra1) + uy*(pair.sigy[i]*fra2 + pair.sigy[i+1]*fra1) + uz*(pair.sigz[i]*fra2 + pair.sigz[i+1]*fra1);

   //add syn signal to syn seismo
   i = 0;
   ne = nb+pair.shdz.npts;
   if(nb<0) { i = -nb; nb = 0; }
   if(ne>npts) ne = npts;
   int iseg, isegb = nb/segnpt, isege = (ne-2)/segnpt;
   for(iseg=isegb;iseg<=isege;iseg++) pthread_mutex_lock(&(lock[ista][iseg]));
   for(;nb<ne-1;i++,nb++) (*sig2)[nb] += sigtmp[i];
   for(iseg=isegb;iseg<=isege;iseg++) pthread_mutex_unlock(&(lock[ista][iseg]));

   //shd.b = 0.;
   //shd.npts = npts;
   *sd = pair.shdz;
   return 1;
}

void ReadSRC(char *fname) {
   FILE *fin;
   char buff[300];
   int i, nblk;
   if ( (fin = fopen(fname, "r")) == NULL ) {
      cout<<"Cannot open file "<<fname<<" to read."<<endl;
      exit(0);
   }
   src = NULL;
   for(i=0,nblk=0 ;fgets(buff, 300, fin)!=NULL; i++) {
      if(i>=nblk*NBLKsrc) src = (struct SRC *)realloc(src, (++nblk)*NBLKsrc*sizeof(struct SRC));
      sscanf(buff, "%f %f %f %f", &(src[i].lon), &(src[i].lat), &(src[i].amp), &(src[i].foccur));
   }
   fclose(fin);
   nsrc = i;
}

void InitSTA(char *fname) {
   FILE *fin;
   char buff[300];
   int i, nblk;
   if ( (fin = fopen(fname, "r")) == NULL ) {
      cout<<"Cannot open file "<<fname<<" to read."<<endl;
      exit(0);
   }
   sta = NULL;
   for(i=0,nblk=0 ;fgets(buff, 300, fin)!=NULL; i++) {
      if(i>=nblk*NBLKsta) sta = (struct STA *)realloc(sta, (++nblk)*NBLKsta*sizeof(struct STA));
      sscanf(buff, "%f %f %s", &(sta[i].lon), &(sta[i].lat), sta[i].name);
      sta[i].sig = (float *) calloc (npts, sizeof(float));
      sta[i].shd.user1 = 0;
   }
   fclose(fin);
   nsta = i;
}

void *PrepareSource(void *rid) {
   char name[100];
   int ista, isrc, nthread;
   int npt;
   float deci;
   double dis, azi;
   char fkpath[100] = "/home/tianye/Software/fk/";

 //  struct iRNG range = *((struct iRNG *)rid);
   nthread = *((int *)rid);

   //for(isrc=range.b;isrc<range.e;isrc++) {
   for(isrc=nthread;isrc<nsrc;isrc+=NTHRDS) {
      cout<<"Thread "<<nthread<<": Computing source-station signals for the "<<isrc+1<<"/"<<nsrc<<"th source"<<endl;
      for(ista=0;ista<nsta;ista++) {
         calc_dist(src[isrc].lat, src[isrc].lon, sta[ista].lat, sta[ista].lon, &dis);
         calc_azimuth(src[isrc].lat, src[isrc].lon, sta[ista].lat, sta[ista].lon, &azi);
         //distance roundoff
         deci = (int)floor(dis/100.)/10.;
         if(deci==0.) dis = floor(dis/0.02+0.5)*0.02;
         else dis = floor(dis/deci+0.5)*deci;
         sprintf(name, "Q100_%d/%.2f.grn.0", dep, dis);
         if( access( name, R_OK) == -1 ) {
            npt = (int) pow( 2., (ceil)(log((dis*0.4+50.)/dt)/log(2.)) );
            sprintf(buff, "%sfk.pl -MQ100/%d/k -N%d/%f %.2f", fkpath, dep, npt, dt, dis);
            system(buff);
         }
         sprintf(name, "SynSig_Q100_%d/dis%.2f_azi%.2f_x.z", dep, dis, azi);
         if( access( name, R_OK) == -1 ) {
            sprintf(buff, "%ssyn -M1/0/0 -D1 -A%f -O%s -GQ100_%d/%.2f.grn.0", fkpath, azi, name, dep, dis);
            system(buff);
         }
         if ( read_sac(name, &(pairs[isrc][ista].sigx), &(pairs[isrc][ista].shdx))==NULL ) cout<<"SAC reading error!"<<endl;
         sprintf(name, "SynSig_Q100_%d/dis%.2f_azi%.2f_y.z", dep, dis, azi);
         if( access( buff, R_OK) == -1 ) {
            sprintf(buff, "%ssyn -M1/90/0 -D1 -A%f -O%s -GQ100_%d/%.2f.grn.0", fkpath, azi, name, dep, dis);
            system(buff);
         }
         if ( read_sac(name, &(pairs[isrc][ista].sigy), &(pairs[isrc][ista].shdy))==NULL ) cout<<"SAC reading error!"<<endl;
         sprintf(name, "SynSig_Q100_%d/dis%.2f_azi%.2f_z.z", dep, dis, azi);
         if( access( buff, R_OK) == -1 ) {
            sprintf(buff, "%ssyn -M1/0/90 -D1 -A%f -O%s -GQ100_%d/%.2f.grn.0", fkpath, azi, name, dep, dis);
            system(buff);
         }
         if ( read_sac(name, &(pairs[isrc][ista].sigz), &(pairs[isrc][ista].shdz))==NULL ) cout<<"SAC reading error!"<<endl;
      }
   }
   pthread_exit(NULL);
}

void *ActSource( void * rid ) {
   int ista, isrc, nthread;
   float tmin = -3000./2., tmax = (npts-1)*dt, tstep = 1.;
   float stk, dip, time;

   nthread = *((int *)rid);

   for(isrc=nthread;isrc<nsrc;isrc+=NTHRDS) {
      if( src[isrc].amp == 0 ) continue;
      cout<<"Thread "<<nthread<<": Working on source "<<isrc+1<<"/"<<nsrc<<endl;
      for(time=tmin; time<tmax; time+=tstep) {
         if( ! RandSRC(src[isrc].foccur*tstep, &stk, &dip) ) continue;
         for(ista=0;ista<nsta;ista++)
            if( ! ClipAddSac(isrc, ista, stk, dip, time, &(sta[ista].sig), &(sta[ista].shd)) ) sta[ista].shd.user1++;
      }
   }
   pthread_exit(NULL);
}


int main(int argc, char *argv[])
{
   if(argc!=3) {
      cout<<"Usage: "<<argv[0]<<" [src infile (lon lat amp(dyne) f_occurrence(Hz))] [station infile (lon lat name)]"<<endl;
      exit(-1);
   }

   //read in noise source locations, amplitudes and frequency of occurrence
   ReadSRC( argv[1]);

   //read in station locations and names
   InitSTA( argv[2]);

   cout<<nsrc<<" sources and "<<nsta<<" stations read in."<<endl;

   //Initialize PAIR struct according to nsrc & nsta
   int i, j;
   pairs = (struct PAIR **) malloc ( nsrc * sizeof(struct PAIR *) );
   for(i=0;i<nsrc;i++) pairs[i] = (struct PAIR *) malloc (nsta * sizeof(struct PAIR));

   /*/divide nsrc sources into NTHRDS groups
   struct iRNG range[NTHRDS];
   int size = (int)ceil((float)nsrc/NTHRDS);
   for(i=0;i<NTHRDS;i++) range[i].b = i*size;
   for(i=0;i<NTHRDS-1;i++) range[i].e = range[i+1].b;
   range[i].e = nsrc; */

   //Thread ids and attributes
   pthread_t tid[NTHRDS];
   pthread_attr_t attr_j;
   pthread_attr_init(&attr_j);
   pthread_attr_setdetachstate(&attr_j, PTHREAD_CREATE_JOINABLE);
   int targs[NTHRDS];

   //initialize mutex(lock) matrix
   lock = (pthread_mutex_t **) malloc ( nsta * sizeof(pthread_mutex_t *) );
   for(i=0;i<nsta;i++) lock[i] = (pthread_mutex_t *) malloc ( segnum * sizeof(pthread_mutex_t) );
   for(i=0;i<nsta;i++) for(j=0;j<segnum;j++) pthread_mutex_init(&(lock[i][j]), NULL);

   //compute/readin synthetic signals between all source-station pairss
   sprintf(buff, "mkdir -p SynSig_Q100_%d", dep);
   system(buff);
   int rc;
   for(i=0;i<NTHRDS;i++) {
      targs[i] = i;
      rc = pthread_create( &tid[i], &attr_j, PrepareSource, (void *) (&(targs[i])) );
      if(rc) {
         cout<<"Creation failed!  ERR: "<<strerror(rc)<<endl;
	 exit(0);
      }
   }
   for(i=0;i<NTHRDS;i++) pthread_join(tid[i], NULL);

   //compute daily summation of signals on each station
   for(i=0;i<NTHRDS;i++) {
      targs[i] = i;
      rc = pthread_create( &tid[i], &attr_j, ActSource, (void *) (&(targs[i])) );
      if(rc) {
         cout<<"Creation failed! ERR: "<<strerror(rc)<<endl;
         exit(0);
      }
   }
   for(i=0;i<NTHRDS;i++) pthread_join(tid[i], NULL);

   //free attributes
   pthread_attr_destroy(&attr_j);

   //output and cleanup for each station
   int ista, isrc;
   for(ista=0;ista<nsta;ista++) {
      sta[ista].shd.npts = npts;
      sta[ista].shd.b = 0.;
      sprintf(buff, "%s.BHZ.SAC", sta[ista].name);
      write_sac(buff, sta[ista].sig, &(sta[ista].shd));
      free(sta[ista].sig);
   }

   //cleanup mutex matrix
   for(i=0;i<nsta;i++) {
      for(j=0;j<segnum;j++) pthread_mutex_destroy(&lock[i][j]);
      free(lock[i]);
   }
   free(lock);

   //cleanup
   for(isrc=0;isrc<nsrc;isrc++) {
      for(ista=0;ista<nsta;ista++) { free(pairs[isrc][ista].sigx); free(pairs[isrc][ista].sigy); free(pairs[isrc][ista].sigz); }
      free(pairs[isrc]);
   }
   free(pairs); free(src); free(sta);

   pthread_exit(NULL);
}
