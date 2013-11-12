#ifndef SACREC_H
#define SACREC_H

#include "SacOps.h"
#include "DisAzi.h"
#include "PathAverage.h"
#include <cmath>


class SACREC {
protected:
   char fsac[150];
   SAC_HD shd;
   float *sig;
   float cper, grv;
public:
   /* initialize by setting sac file name, center period and group speed */
   SACREC( const char *inname = NULL, const float cperin = 0., const float grvin = 0. ) {
      sig = NULL;
      sprintf(fsac, "%s", inname);
      cper = cperin; grv = grvin;
   }
   void setGrv( float &grvin ) { grv = grvin; }

   /* load and access sac header and signal */
   bool load() { 
      if( ! read_sac(fsac, &sig, &shd) ) return false; 
      if( shd.depmin != shd.depmin || shd.depmax != shd.depmax ) return false;
      if( shd.dist <= 0. ) { calc_dist(shd.evla, shd.evlo, shd.stla, shd.stlo, &shd.dist); }
      return true;
   }

   const SAC_HD &GetHeader() { return shd; }

   /* compute aamplitude and rms noise */
   float amp(float time) { 
      if( sig == NULL ) return -12345.; 
      int i = (int)floor((time-shd.b)/shd.delta+0.5);
      if( i < 0 || i >= shd.npts ) return -12345.;
      return sig[i]; 
   }

   float noise() {
      if( sig == NULL ) { return -12345.; }
      float maxT, minT = shd.dist/grv + cper*5. + 500.; //move out of the signal tail
      if( (shd.e - minT) < 50. ) { // no enough noise for computing rms
         return -1.;
      }
      else if( (shd.e - minT) < 600. ) { maxT = shd.e - 10.; }
      else {
         minT = shd.e - 600.;
         maxT = shd.e - 100.;
      }
      // compute rms noise
      int ib = (int)floor(minT/shd.delta);
      int ie = (int)ceil(maxT/shd.delta);
      float noiserms=0.;
      for(int i=ib;i<ie;i++) noiserms += sig[i] * sig[i];
      noiserms=sqrt(noiserms/(ie-ib-1.));
      return noiserms;
   }

   /* station names, locations, and distance */
   char * Sta1() { return shd.kevnm; }
   char * Sta2() { return shd.kstnm; }
   Point<float> P1() { return Point<float>(shd.evlo, shd.evla); }
   Point<float> P2() { return Point<float>(shd.stlo, shd.stla); }
   double Dist() { return shd.dist; }

   /* search for maximum amplitude of precursoring signal */
   bool Precursor( float *amp, float *time ) {
      if( sig == NULL ) { return false; }
      float Tmax = shd.dist/grv - cper;
      int ib = (int)floor((-Tmax - shd.b)/shd.delta + 1.5);
      if( ib < 1 ) ib = 1;
      int ie = (int)floor((Tmax - shd.b)/shd.delta + 0.5);
      if( ie > shd.npts-1 ) ie = shd.npts-1;
      /* check for largest local maximum */
      float sigtmp[3], precamp = 0., prectime;
      for(int i=ib; i<ie; i++) {
	 for(int j=0; j<3; j++) { sigtmp[j] = sig[i-1+j]; if(sigtmp[j]<0.) sigtmp[j] = -sigtmp[j];}
	 if( (sigtmp[1]-sigtmp[0]) <= 0. || (sigtmp[1]-sigtmp[2]) <= 0. ) continue;
	 if( sigtmp[1] > precamp ) { precamp = sigtmp[1]; prectime = i; }
      }
      prectime = shd.b + prectime*shd.delta;
      *amp = precamp; *time = prectime;
      if( precamp > 0. ) return true;
      else return false;
   }

   ~SACREC() { if( sig ) delete [] sig; }
};

#endif
