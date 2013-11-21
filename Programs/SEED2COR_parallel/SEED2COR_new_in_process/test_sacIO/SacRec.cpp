#include "SacRec.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
//#include <pthread.h>

/*
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

pthread_mutex_t fiolock;
int jday ( int y, int m, int d ) {
   int i, jd = 0;
   for( i = 1; i < m; i++ ) {
      if ( (i==1) || (i==3) || (i==5) || (i==7) || (i==8) || (i==10) ) jd += 31;
      else if (i==2) {
         if ( (y%400==0) || (y%100!=0&&y%4==0) ) jd += 29;
         else jd += 28;
      }
      else jd += 30;
   }
   return jd + d;
}

double abs_time ( int yy, int jday, int hh, int mm, int ss, int ms ) {
     //computes time in s relative to 1900
   int nyday = 0, i;
   for( i = 1901; i < yy; i++ ) {
      if ( (i%400==0) || (i%100!=0&&i%4==0) ) nyday += 366;
      else nyday += 365;
   }
   return 24.*3600.*(nyday+jday) + 3600.*hh + 60.*mm + ss + 0.001*ms;
}
*/

/* reformat header time if shd.nzmsec is modified and is out of the range [0,1000) */
void SacRec::UpdateTime() {
   if(shd.nzmsec < 1000 && shd.nzmsec>=0) return;
   int i = (int)floor(shd.nzmsec/1000);
   shd.nzmsec -= i*1000;
   shd.nzsec += i;
   if(shd.nzsec < 60 && shd.nzsec>=0) return;
   i = (int)floor(shd.nzsec/60);
   shd.nzsec -= 60*i;
   shd.nzmin += i;
   if(shd.nzmin < 60 && shd.nzmin>=0) return;
   i = (int)floor(shd.nzmin/60);
   shd.nzmin -= i*60;
   shd.nzhour += i;
   if(shd.nzhour < 24) return;
   shd.nzhour -= 24;
   shd.nzjday++;
   if(shd.nzjday < 366) return;
   if( ((shd.nzyear%400==0) || (shd.nzyear%100!=0 && shd.nzyear%4==0)) && shd.nzjday<367 ) return;
   shd.nzjday = 1;
   shd.nzyear++;
}

/* load sac header from file 'fname' */
bool SacRec::LoadHD () {
   std::ifstream fsac(fname.c_str());
   if( ! fsac ) return false;
   //pthread_mutex_lock(&fiolock);
   fsac.read( reinterpret_cast<char *>(&shd), sizeof(SAC_HD) );
   fsac.close();
   //pthread_mutex_unlock(&fiolock);
   
   return true;
}

/* load sac header+signal from file 'fname', memory is allocated on heap */
bool SacRec::Load () {
   std::ifstream fsac(fname.c_str());
   if( ! fsac ) return false;
   //pthread_mutex_lock(&fiolock);
   fsac.read( reinterpret_cast<char *>(&shd), sizeof(SAC_HD) );
   //sig = std::make_shared<float>( new float[shd.npts*sizeof(float)] );
   float* sigtmp = new float[shd.npts];
   fsac.read( reinterpret_cast<char *>(sigtmp), sizeof(float)*shd.npts );
   fsac.close();
   sig = std::unique_ptr<float[]>(sigtmp);
   //pthread_mutex_unlock(&fiolock);

   /* calculate t0 */
   char koo[9];
   for ( int i=0; i<8; i++ ) koo[i] = shd.ko[i]; koo[8] = 0;
   float fes;
   int eh, em;
   sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);
   shd.o = shd.b + (shd.nzhour-eh)*3600. + (shd.nzmin-em)*60. + shd.nzsec-fes + shd.nzmsec*.001;

   return true;
}

/* write to file '*outfname' */
bool SacRec::Write (const char *outfname) {
   /* check if signal is loaded */
   if( ! sig ) {
      std::cerr<<"ERROR(write_sac): No signal loaded in the memory! "<<outfname<<std::endl;
      return false;
   }
   /* open file */
   std::ofstream fsac(outfname);
   if( ! fsac ) {
      std::cerr<<"ERROR(write_sac): Cannot open file "<<outfname<<std::endl;
      return false;
   }
   /* update header */
   shd.iftype = (int)ITIME;
   shd.leven = (int)TRUE;
   shd.lovrok = (int)TRUE;
   shd.internal4 = 6L;

  /* search min and max amplitude */
   shd.depmin = sig[0];
   shd.depmax = sig[0];
   for ( int i = 0; i < shd.npts ; i++ ) {
       if ( shd.depmin > sig[i] ) shd.depmin = sig[i];
       else if ( shd.depmax < sig[i] ) shd.depmax = sig[i];
   }

   /* check and re-format header time if necessary */
   UpdateTime();

   //pthread_mutex_lock(&fiolock);
   fsac.write( reinterpret_cast<char *>(&shd), sizeof(SAC_HD) );
   fsac.write( reinterpret_cast<char *>(sig.get()), sizeof(float)*shd.npts );
   //pthread_mutex_unlock(&fiolock);

   fsac.close();
}

bool SacRec::MinMax () {

}
/*
int read_rec(int rec_flag, char *fname, int len, int *rec_b, int *rec_e, int *nrec) {
   FILE *frec;
   int irec;
   if( rec_flag ) {
      if((frec = fopen(fname,"r")) == NULL) return 0;
      pthread_mutex_lock(&fiolock);
      for(irec=0;;irec++)
         if(fscanf(frec,"%d %d", &rec_b[irec], &rec_e[irec])!=2) break;
      *nrec=irec;
      fclose(frec);
      pthread_mutex_unlock(&fiolock);
      if(irec==0) return 0;
   }
   else {
      rec_b[0]=0; rec_e[0]=len-1;
      *nrec=1;
   }
   return 1;
}
*/
