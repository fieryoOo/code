#ifndef SACREC_H
#define SACREC_H

#include "mysac64.h"
#include <cstdio>
#include <string>
#include <iostream>
#include <memory>

#ifndef HAVE_NULLPTR
/* a workaround for nullptr when compiled by old gcc compiler */
const                        // this is a const object...
class {
public:
  template<class T>          // convertible to any type
    operator T*() const      // of null non-member
    { return 0; }            // pointer...
  template<class C, class T> // or any type of null
    operator T C::*() const  // member pointer...
    { return 0; }
private:
  void operator&() const;    // whose address can't be taken
} nullptr = {};              // and whose name is nullptr
#endif

class SacRec {
public:
   std::string fname;
   SAC_HD shd;
   std::unique_ptr<float[]> sig;
   //std::auto_ptr<float> sig;
public:
   /* constructors */
   SacRec( const char* fnamein) : fname(fnamein), sig(nullptr), shd(sac_null) {}
   SacRec( const SacRec& recin )
    : fname(recin.fname), shd(recin.shd), sig(new float[recin.shd.npts]) { 
      std::copy(recin.sig.get(), recin.sig.get()+recin.shd.npts, sig.get()); 
   }
   /* operators */
   SacRec &operator= ( const SacRec& recin ) { 
      fname = recin.fname; shd = recin.shd;
      int npts=recin.shd.npts; sig.reset(new float[npts]); 
      std::copy(recin.sig.get(), recin.sig.get()+npts, sig.get()); 
   }
   /* load sac header from file 'fname' */
   bool LoadHD ();
   /* read sac header+signal from file 'fname', memory is allocated on heap */
   bool Load ();
   /* write to file '*fname' */
   bool Write ( const char *fname );
   /* update/reformat header time if shd.nzmsec is modified and is out of the range [0,1000) */
   void UpdateTime();
   /* search for min and max signal amplitudes */
   bool MinMax();
};


#endif
