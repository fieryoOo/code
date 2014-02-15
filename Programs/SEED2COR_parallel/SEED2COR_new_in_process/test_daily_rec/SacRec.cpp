#include "SacRec.h"
#include <fftw3.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <limits>
#include <chrono>
#include <random>
//#include <pthread.h>


/* ---------------------------------------- Pimpl handle struct ---------------------------------------- */
struct SacRec::SRimpl {

   /* ---------- FFT operations ---------- */
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

   float FillGap( float *pbeg, float *pend, float amp, int hlen, int step ) {
      // defube random number generator
      unsigned timeseed = std::chrono::system_clock::now().time_since_epoch().count();
      std::default_random_engine generator (timeseed);
      std::uniform_real_distribution<float> distribution(-amp, amp);
      auto rand = std::bind ( distribution, generator );
      // parameters
      float oostep = 1./step, slope;
      float *p, *pmid = pbeg + reinterpret_cast<long>( (pend-pbeg) / 2 );
      float alpha = -0.5 / (hlen * hlen);
      // generate tapered random numbers for the 1st half
      *pbeg = rand();
      for(p=pbeg+step; p<pmid; p+=step) { 
	 int ndiff = (p-pbeg);
	 float gdamp = exp( alpha * ndiff * ndiff );
	 *p = rand() * gdamp;
	 slope = ( *p - *(p-step) ) * oostep;
	 for(int i=1; i<step; i++) *(p+i) = *p + slope * i;
      }
      // generate tapered random numbers for the 2nd half
      for(; p<pend; p+=step) { 
	 int ndiff = (pend-p);
	 float gdamp = exp( alpha * ndiff * ndiff );
	 *p = rand() * gdamp; 
	 slope = ( *p - *(p-step) ) * oostep;
	 for(int i=1; i<step; i++) *(p+i) = *p + slope * i;
      }
      // connect the last several points
      *pend = rand();
      p = p-step;
      slope = ( *pend - *p ) / ( pend - p );
      for(p=p+1; p<pend; p++) *p = *(p-1) + slope;

   }

   float av_sig (float *sig, int i, int N, int nwin ) {
      int n1, n2, j, nav = 0;
      float av = 0.;

      if ( nwin > N ) nwin = N;
      n1 = i - nwin/2;
      if ( n1 < 0 ) n1 = 0;
      n2 = n1 + nwin - 1;
      if ( n2 > N-1 ) n2 = N-1;
      n1 = n2 - nwin + 1;

      for ( j = n1; j <= n2; j++ ) 
         if ( sig[j] < 1.e29 ) {
            av += sig[j];
            nav++;
         }

      if ( nav < 1 ) av = std::numeric_limits<float>::max();
      else av = av/(float)nav;

      return av;
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

/* ---------- sac operations ---------- */
//#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
//#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

//pthread_mutex_t fiolock;

double SacRec::AbsTime () {
   if( ! sig ) return -1.;
   //computes time in s relative to 1900
   int nyday = 0;
   for( int i=1901; i<shd.nzyear; i++ ) {
      if ( (i%400==0) || (i%100!=0&&i%4==0) ) nyday += 366;
      else nyday += 365;
   }
   return 24.*3600.*(nyday+shd.nzjday) + 3600.*shd.nzhour + 60.*shd.nzmin + shd.nzsec + 0.001*shd.nzmsec;
}


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
bool SacRec::RMSAvg ( float tbegin, float tend, int step, float& rms ) {
   if( ! sig ) return false;

   rms = 0.;
   float maxfloat = std::numeric_limits<float>::max();
   int neff = 0, ibeg = nint((tbegin-shd.b)/shd.delta), iend = nint((tend-shd.b)/shd.delta+1.);
   for ( int i = ibeg; i < iend ; i+=step ) {
      if( sig[i] >= maxfloat ) continue;
      rms += sig[i] * sig[i];
      neff++;
   }
   rms = sqrt(rms/(neff-1));
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


/* ---------------------------------------- cut and merge ---------------------------------------- */

bool SacRec::merge( SacRec sacrec2 ) {
   // make sure that both signals are loaded
   if( !sig || !sacrec2.sig ) return false;
   
   SAC_HD& shd2 = sacrec2.shd;
   // starting and ending time
   double t1b = this->AbsTime();
   double t2b = sacrec2.AbsTime();
   double t1e = t1b + (shd.npts-1)*shd.delta;
   double t2e = t2b + (shd2.npts-1)*shd2.delta;
   double T1 = std::min(t1b, t2b), T2 = std::max(t1e, t2e);

   double dt = (int)floor(shd.delta*1e8+0.5)/1e8, tshift;
   int N = (int)floor((T2-T1)/shd.delta+0.5)+1;

   /* pre-merging checks */
   if( N > 2.6e7 ) {
      std::cerr<<"Error(SacRec::merge): The merged signal will be larger than 100M. Stopped! " <<std::endl;
      return false;
   }
   if( (shd.delta-shd2.delta) > 0.0001 ) {
      std::cerr<<"Error(SacRec::merge): sps mismatch! Stopped! " <<std::endl;
      return false;
   }

   /* allocate new space */
   std::unique_ptr<float[]> sig0(new float[N]);
   std::fill(&(sig0[0]), &(sig0[N]), std::numeric_limits<float>::max()); // initialize the array to max float
   //for (j=0;j<N;j++) sig0[j] = 1.e30;
   //std::copy(&(sig1[0]), &(sig1[0])+5, &(sig0[0]));

   /* compute merge locations */
   int nb;
   bool reversed;
   if( t1b > t2b ) {
      reversed = true;
      nb = (int)floor((t1b-t2b)/dt+0.5);
      tshift = (shd.b-shd2.b) + (nb*dt-(t1b-t2b));
   }
   else {
      reversed = false;
      nb = (int)floor((t2b-t1b)/dt+0.5);
      tshift = (shd2.b-shd.b) + (nb*dt-(t2b-t1b));
   }
   if( fabs(tshift) > 1.e-3 )
      std::cerr << "Warning(SacRec::merge): signal shifted by " << tshift << "sec when merging ("<< fname.c_str() << ") !" << std::endl;


   /* copy in the signals */
   if( reversed ) {
      std::copy(&(sacrec2.sig[0]), &(sacrec2.sig[0])+shd2.npts, &(sig0[0]));
      std::copy(&(sig[0]), &(sig[0])+shd.npts, &(sig0[nb]));
   }
   else {
      std::copy(&(sig[0]), &(sig[0])+shd.npts, &(sig0[0]));
      std::copy(&(sacrec2.sig[0]), &(sacrec2.sig[0])+shd2.npts, &(sig0[nb]));
   } 

   /* resign sig pointer and update sac header */
   sig = std::move(sig0);
   if( reversed ) shd = shd2;
   shd.npts = N;
}

int SacRec::arrange(const char* recname) {
   // count holes
   float maxfloat = std::numeric_limits<float>::max();
   int Nholes=0;
   for(int i=0;i<shd.npts;i++)
      if(sig[i]>=maxfloat) Nholes++;

   // produce record file for holes if a recname is given
   if( recname ) {
      int rec_b[1000], rec_e[1000];
      //char recname[200];
      rec_b[0] = 0;
      int j=0;
      for(int i=1;i<shd.npts;i++) {
         if(sig[i-1] >= maxfloat) { if(sig[i] < maxfloat) rec_b[j] = i; }
         else if(sig[i] >= maxfloat) rec_e[j++] = i;
      }
      if(sig[shd.npts-1]<maxfloat) rec_e[j++] = shd.npts;
      //sprintf(recname, "%s_rec1", ffsac.c_str());
      std::ofstream frec(recname);
      for(int i=0;i<j;i++) frec << rec_b[i] << " " << rec_e[i] << std::endl;
      //fprintf(frec, "%d %d\n", rec_b[i], rec_e[i]);
      frec.close();
   }

   // fill gaps with random numbers
   int ib, hlen = 100./shd.delta, step = std::min(0.5/shd.delta, 1.);
   float sigrms;
   this->RMSAvg(shd.b, shd.e, 50./shd.delta, sigrms);
   bool isgap = false;
   for (int i=0; i<shd.npts; i++ ) {
      if( isgap ) {
	 if ( sig[i] < maxfloat ) {
	    pimpl->FillGap(&(sig[ib]), &(sig[i]), sigrms, hlen, step );
	    isgap = false;
	 }
      }
      else {
	 if ( sig[i] >= maxfloat ) {
	    ib = i;
	    isgap = true;
	 }
      }
   }
   if( isgap ) pimpl->FillGap(&(sig[ib]), &(sig[0])+shd.npts, sigrms, hlen, step );

   /*
   float av;
   int npart;
   for (int i=0; i<shd.npts; i++ ) if ( sig[i] > 1.e29 ) {
      for(npart=256;npart>1;npart*=0.5) {
         av = pimpl->av_sig (&(sig[0]), i, shd.npts, shd.npts/npart );
         if ( av < 1.e29 ) break;
      }
      if(npart==1) av=0.;
      sig[i] = av;
   }
   */

   return Nholes;
}


bool SacRec::Resample( float sps ) {
   if( ! sig ) return false;
   int ithread = 0;
   /* anti-aliasing filter */
   float dt = 1./sps;
   int iinc = (int)floor(dt/shd.delta+0.5);
   if(iinc!=1) {
      double f1 = -1., f2 = 0., f3 = sps/2.2, f4 = sps/2.01;
      Filter(-1., 0., f3, f4);
   }

   /* allocate space for the new sig pointer */
   int i, j;
   int nptst = (int)floor((shd.npts-1)*shd.delta*sps+0.5)+10;
   float nb;
   std::unique_ptr<float[]> sig2(new float[nptst]);
   //if( (*sig2 = (float *) malloc (nptst * sizeof(float)))==NULL ) perror("malloc sig2");
   long double fra1, fra2;
   nb = ceil((shd.nzmsec*0.001+shd.b)*sps);
   i = (int)floor((nb*dt-shd.nzmsec*0.001-shd.b)/shd.delta);
   if(fabs(iinc*shd.delta-dt)<1.e-7) { //sps is a factor of 1/delta
      fra2 = (nb*dt-i*shd.delta-shd.nzmsec*0.001-shd.b)/shd.delta;
      fra1 = 1.-fra2;
      if(fra2==0)
         for(j=0;i<shd.npts;j++) {
            sig2[j] = sig[i];
            i += iinc;
         }
      else
         for(j=0;i<shd.npts-1;j++) {
            sig2[j] = sig[i]*fra1 + sig[i+1]*fra2;
            i += iinc;
         }
   }
   else { //sps isn't a factor, slower way
      //reports << "*** Warning: sps isn't a factor of " << (int)floor(1/shd.delta+0.5) << ", watch out for rounding error! *** ";
      std::cerr << "Warning(SacRec::Resample): sps isn't a factor of " << (int)floor(1/shd.delta+0.5) << ", watch out for rounding error ("<<fname.c_str()<<") !"<<std::endl;
      long double ti, tj;
      iinc = (int)floor(dt/shd.delta);
      ti = i*shd.delta+shd.nzmsec*0.001+shd.b;
      tj = nb*dt;
      for(j=0;i<shd.npts-1;j++) {
         fra2 = tj-ti;
         sig2[j] = sig[i] + (sig[i+1]-sig[i])*fra2;
         tj += dt;
         i += iinc;
         ti += iinc*shd.delta;
         if( ti+shd.delta <= tj ) { ti += shd.delta; i++; }//if(j%1000==0)cerr<<i<<" "<<ti<<j<<" "<<tj<<" "<<endl;}
      }
   }
   sig = std::move(sig2);
   shd.nzmsec = (int)(nb*dt*1000+0.5);
   shd.b = 0.;
   if(shd.nzmsec>=1000) UpdateTime();
   shd.delta = dt;
   shd.npts = j;
   return true;
}

