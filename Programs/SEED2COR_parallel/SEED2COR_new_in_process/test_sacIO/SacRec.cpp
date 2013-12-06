#include "SacRec.h"
#include <fftw3.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
//#include <pthread.h>


/* ---------------------------------------- Pimpl handle struct ---------------------------------------- */
struct SacRec::SRimpl {
   #define PI 3.14159265358979323846
   // forward FFT 
   void FFTW_F(fftw_plan plan, fftw_complex *out, int ns, float *seis, int n) {
      fftw_execute(plan);
      //pthread_mutex_lock(&fftlock);
      fftw_destroy_plan(plan);
      //pthread_mutex_unlock(&fftlock);
      int k;
      for(k=0; k<n; k++) seis[k] = out[k][0];
      //for(k=0; k<n; k+=1000) if(seis[k] != 0.) printf("%d %f\n", k, seis[k]);
   }
   // backward FFT
   void FFTW_B(int type, float *seis, int n, fftw_complex **in, fftw_complex **out, int *nso, fftw_plan *planF, int Fflag) {
      int ns = (int)(log((double)n)/log(2.))+1;
      if(ns<13) ns = 13;
      ns = (int)pow(2,ns); *nso = ns;
      *in = (fftw_complex *) fftw_malloc ( ns * sizeof(fftw_complex) );//fftw_alloc_complex(ns);
      *out = (fftw_complex *) fftw_malloc ( ns * sizeof(fftw_complex) );//fftw_alloc_complex(ns);

      //measure plan using the allocated in/out blocks
      //pthread_mutex_lock(&fftlock);
      fftw_plan plan = fftw_plan_dft_1d (ns, *in, *out, FFTW_BACKWARD, type); //FFTW_ESTIMATE / FFTW_MEASURE
      if( Fflag == 1 ) *planF = fftw_plan_dft_1d (ns, *out, *in, FFTW_FORWARD, type);
      if( plan==NULL || (Fflag==1 && *planF==NULL) ) {
         fprintf(stderr,"Error(FFTW_B): fftw_plan creation failed!!\n");
         //pthread_mutex_unlock(&fftlock); exit(0);
      }
      //pthread_mutex_unlock(&fftlock);
      //initialize input array and excute
      memset(*in, 0, ns*sizeof(fftw_complex));
      int k;
      for(k=1; k<n; k++) (*in)[k][0] = seis[k];
      fftw_execute(plan);
      //cleanup
      //pthread_mutex_lock(&fftlock);
      fftw_destroy_plan(plan);
      //pthread_mutex_unlock(&fftlock);
      //if( Fflag==0 ) fftw_free(*in);

      //kill half spectrum and correct ends
      int nk = ns/2+1;
      for(k=nk;k<ns;k++) {
         (*out)[k][0] = 0.;
         (*out)[k][1] = 0.;
      }
      (*out)[0][0] /= 2.; (*out)[0][1] /= 2.;
      (*out)[nk-1][1] = 0.;
   }
   // lowpass taper
   void TaperL( double f3, double f4, double dom, int nk, fftw_complex *sf ) {
      double f, ss;
   
      int i = (int)ceil(f3/dom);
      for(f=i*dom; f<f4; i++, f+=dom) {
         ss = ( 1. + cos(PI*(f3-f)/(f4-f3)) ) / 2.;
         sf[i][0] *= ss;
         sf[i][1] *= ss;
      }
      for(;i<nk;i++) {
         sf[i][0] = 0.;
         sf[i][1] = 0.;
      }

      return;
   }
   // bandpass taper
   void TaperB( double f1, double f2, double f3, double f4, double dom, int nk, fftw_complex *sf ) {
      int i;
      double f, ss;

      for(i=0, f=0.; f<f1; i++, f+=dom) {
         sf[i][0] = 0.;
         sf[i][1] = 0.;
      }
      for(; f<f2; i++, f+=dom) {
         ss = ( 1. - cos(PI*(f1-f)/(f2-f1)) ) / 2.;
         sf[i][0] *= ss;
         sf[i][1] *= ss;
      }
      i = (int)ceil(f3/dom);
      for(f=i*dom; f<f4; i++, f+=dom) {
         ss = ( 1. + cos(PI*(f3-f)/(f4-f3)) ) / 2.;
         sf[i][0] *= ss;
         sf[i][1] *= ss;
      }
      for(;i<nk;i++) {
         sf[i][0] = 0.;
         sf[i][1] = 0.;
      }

      return;
   }
   // gaussian taper
   void TaperGaussian( double fcenter, double fhlen, double dom, int nk, fftw_complex *sf ) {
      int i;
      float gauamp;
      double f, fstart, fend, fmax = (nk-1)*dom;
      // define effective window for given gaussian halflength
      f = fhlen * 4.;
      fstart = fcenter - f;
      fend = fcenter + f; 
      if( fend > fmax ) fend = fmax;
      //std::cerr<<"fcenter "<<fcenter<<"  dom "<<dom<<"  fstart "<<fstart<<"  fend "<<fend<<"  fmax "<<fmax<<std::endl;
      // cut high frequecies
      if( fstart > 0. ) {
         for(i=0, f=0.; f<fstart; i++, f+=dom) {
	    sf[i][0] = 0.;
	    sf[i][1] = 0.;
      }
	 f = i * dom; // correct for round-off error
      }
      else { f = 0.; i = 0; }
      // apply taper
      float alpha = -0.5/(fhlen*fhlen);
      for(; f<fend-1.e-10; i++, f+=dom) {
	 gauamp = f - fcenter;
	 gauamp = exp( alpha * gauamp * gauamp );
	 sf[i][0] *= gauamp;
	 sf[i][1] *= gauamp;
      }
      // cut low frequencies
      if( fend < fmax ) {
	 f = i * dom; // again, correct for round-off
	 for(; i<nk; i++) {
	    sf[i][0] = 0.;
	    sf[i][1] = 0.;
	 }
      }
   }

};


/* ---------------------------------------- constructors and operators ---------------------------------------- */
/* default constructor */
SacRec::SacRec( const char* fnamein )
 : sig(nullptr), shd(sac_null), pimpl(new SRimpl() ) {
   if( fnamein ) fname = fnamein;
}

/* copy constructor */
SacRec::SacRec( const SacRec& recin )
 : fname(recin.fname), shd(recin.shd), sig(new float[recin.shd.npts]), pimpl( new SRimpl(*(recin.pimpl)) ) { 
   std::copy(recin.sig.get(), recin.sig.get()+recin.shd.npts, sig.get()); 
}

/* assignment operator */
SacRec& SacRec::operator= ( const SacRec& recin ) { 
   pimpl.reset( new SRimpl(*(recin.pimpl)) );
   fname = recin.fname; shd = recin.shd;
   int npts=recin.shd.npts; sig.reset(new float[npts]); 
   std::copy(recin.sig.get(), recin.sig.get()+npts, sig.get()); 
}

/* destructor */
SacRec::~SacRec() {}

/* ---------------------------------------- sac IO ---------------------------------------- */
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


/* ---------------------------------------- sac operations ---------------------------------------- */
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


/* search for min&max signal positions and amplitudes */
inline static int nint( float in ) { return static_cast<int>(floor(in+0.5)); }
bool SacRec::MinMax (float tbegin, float tend, float& tmin, float& min, float& tmax, float& max) {
   if( ! sig ) return false;
   min = max = sig[0];
   int imin = 0, imax = 0;
   for ( int i = nint((tbegin-shd.b)/shd.delta); i < nint((tend-shd.b)/shd.delta+1.) ; i++ ) {
       if ( min > sig[i] ) { min = sig[i]; imin = i; }
       else if ( max < sig[i] ) { max = sig[i]; imax = i; }
   }
   tmin = shd.b + imin*shd.delta;
   tmax = shd.b + imax*shd.delta;
   return true;
}


/* compute the root-mean-square average in a given window */
bool SacRec::RMSAvg ( float tbegin, float tend, float& rms ) {
   if( ! sig ) return false;
   int i, ibeg = nint((tbegin-shd.b)/shd.delta);
   for ( rms=0.,i = ibeg; i < nint((tend-shd.b)/shd.delta+1.) ; i++ )
      rms += sig[i] * sig[i];
   rms = sqrt(rms/(i-ibeg-1.));
   return true;
}


/* ---------------------------------------- time - frequency ---------------------------------------- */
/* 3 different types of filters */
bool SacRec::Filter ( double f1, double f2, double f3, double f4, SacRec& srout ) {
   if( ! sig ) return false;	// check signal
   if( &srout != this ) {	// initialize srout if not filtering in place 
      srout.fname = fname; srout.shd = shd; 
      srout.sig.reset( new float[shd.npts] );
      srout.pimpl.reset( new SRimpl(*(pimpl)) );
   }
   double dt = shd.delta;
   if(f4 > 0.5/dt) {
      std::cerr<<"Warning(SacRec::Filter): filter band out of range!"<<std::endl;;
      f4 = 0.49999/dt;
   }
   fftw_plan planF = NULL;

   //backward FFT: s ==> sf
   int ns, n=shd.npts;
   fftw_complex *s, *sf;
   pimpl->FFTW_B(FFTW_ESTIMATE, &(sig[0]), n, &s, &sf, &ns, &planF, 1);

   //make tapering
   int nk = ns/2+1;
   double dom = 1./dt/ns;
   if( (f1==-1. || f2==-1.) && (f3>0. && f4>0.) ) pimpl->TaperL( f3, f4, dom, nk, sf );
   else if( f1>=0. && f2>0. && f3>0. && f4>0. ) pimpl->TaperB( f1, f2, f3, f4, dom, nk, sf );
   else if( f1==-1. && f4==-1. ) pimpl->TaperGaussian( f2, f3, dom, nk, sf );
   else {
      std::cerr<<"Warning(SacRec::Filter): Unknow filter type for parameters = "<<f1<<" "<<f2<<" "<<f3<<" "<<f4<<std::endl;
      return false;
   }

   //forward FFT: sf ==> s
   pimpl->FFTW_F(planF, s, ns, &(srout.sig[0]), n);
   fftw_free(s); fftw_free(sf);

   //forming final result
   float ftmp = 2./ns;
   for(int k=0; k<n; k++) srout.sig[k] *= ftmp;
   //   if( seis_in[k]==0 ) seis_out[k] = 0.;
   //   else seis_out[k] *= ftmp;

   return true;
}


