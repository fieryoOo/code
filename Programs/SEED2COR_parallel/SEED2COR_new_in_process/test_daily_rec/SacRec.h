#ifndef SACREC_H
#define SACREC_H

#include "mysac64.h"
#include <cstdio>
#include <cstddef>
#include <string>
#include <iostream>
#include <memory>
/*
// a workaround for nullptr when compiled by old gcc compiler
const class {				// this is a const object...
public:
   template<class T>			// convertible to any type
   operator T*() const { return 0 }	// of null non-member pointer...
   template<class C, class T>		// or any type of null
   operator T C::*() const { return 0; }// member pointer...
private:
  void operator&() const;		// whose address can't be taken
} nullptr = {};				// and whose name is nullptr
*/


class SacRec {
private:
   /* impl pointer */
   struct SRimpl;
   std::unique_ptr<SRimpl> pimpl;
public:
   std::string fname;			// input file name
   SAC_HD shd;				// sac header
   std::unique_ptr<float[]> sig;	// pointer the the signal
   //std::auto_ptr<float> sig;
public:
   /* constructors */
   SacRec( const char* fnamein = NULL );	// default
   SacRec( const SacRec& recin );		// copy
   /* operators */
   SacRec &operator= ( const SacRec& recin );	// assignment
   /* load sac header from file 'fname' */
   bool LoadHD( const char* fnamein ) { if( fnamein ) fname = fnamein; return LoadHD(); }
   bool LoadHD ();
   /* read sac header+signal from file 'fname', memory is allocated on heap */
   bool Load( const char* fnamein ) { if( fnamein ) fname = fnamein; return Load(); }
   bool Load ();
   /* write to file '*fname' */
   bool Write ( const char *fname );
   /* compute the absolute time in sec relative to 1900.01.00 */
   double AbsTime();
   /* update/reformat header time if shd.nzmsec is modified and is out of the range [0,1000) */
   void UpdateTime();
   /* search for min&max signal positions and amplitudes */
   bool MinMax ( float tbegin, float tend, float& tmin, float& min, float& tmax, float& max );
   /* compute the root-mean-square average in a given window */
   bool RMSAvg ( float tbegin, float tend, float& rms) { return RMSAvg( tbegin, tend, 1, rms); }
   bool RMSAvg ( float tbegin, float tend, int step, float& rms);
   /* merge a second sacrec to the current */
   bool Merge( SacRec sacrec2 ) {
      merge( sacrec2 );
      arrange();
   }
   bool merge( SacRec sacrec2 );
   int arrange( const char *recname = NULL );
   /* method that performs 3 different types of filters
    * lowpass when ( (f1==-1. || f2==-1.) && (f3>0. && f4>0.) )
    * bandpass when ( f1>=0. && f2>0. && f3>0. && f4>0. )
    * gaussian when ( f1==-1. && f4==-1. ) where f2 = center freqency and f3 = frequency half length */
   bool Filter ( double f1, double f2, double f3, double f4 ) { return Filter(f1, f2, f3, f4, *this); }	// in-place
   bool Filter ( double f1, double f2, double f3, double f4, SacRec& srout );				// out-of-place
   /* resample (with anti-aliasing filter) the signal to given sps */
   bool Resample( float sps );
   /* destructor */
   ~SacRec(); 
};

/*
class AMPH {
private:
   class APimpl;
   std::unique_ptr<APimpl> pimpl;
   std::unique_ptr<float[]> am, ph;
};
*/

#endif
