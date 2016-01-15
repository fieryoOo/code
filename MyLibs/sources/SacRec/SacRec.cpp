#include "SacRec.h"
#include "DisAzi.h"
#include "Parabola.h"
//#include "MyLogger.h"
//#include "SysTools.h"
#include <fftw3.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <limits>
#include <chrono>
#include <random>
#include <algorithm>
//#include <pthread.h>

//#include "SysTools.h"
//extern MEMO memo;

/* ---------------------------------------- Pimpl handle struct ---------------------------------------- */
struct SacRec::SRimpl {

	/* ---------- search by filename in ${PATH} ---------- */
   bool FindInPath( const std::string fname, std::string& absname ) {
      char* pPath = getenv("PATH");
      if( ! pPath ) return false;
      std::stringstream sPath(pPath);
      for(std::string pathcur; std::getline(sPath, pathcur, ':'); ) {
         std::string testname = pathcur + '/' + fname;
         std::ifstream fin(testname);
         if( fin ) {
            absname = testname;
            return true;
         }
      }
      return false;
   }

   /* ---------- FFT operations ---------- */
   //#define PI 3.14159265358979323846
   // forward FFT: out ==> seis
   void FFTW_F(fftw_plan plan, fftw_complex *out, float *seis, int n, const short outtype = 0) {
		/* outtype: real(0)/imaginary(1)/amp(2)/phase(3) of the IFFT result */
      fftw_execute(plan);
      //pthread_mutex_lock(&fftlock);
		#pragma omp critical(fftw)
		{
      fftw_destroy_plan(plan);
		}
      //pthread_mutex_unlock(&fftlock);
      int k;
		switch( outtype ) {
			case 1:
				for(k=0; k<n; k++) seis[k] = out[k][1];
				break;
			case 2:
				for(k=0; k<n; k++) seis[k] = hypot(out[k][0], out[k][1]);	//sqrt( out[k][0]*out[k][0] + out[k][1]*out[k][1] );
				break;
			case 3:
				for(k=0; k<n; k++) seis[k] = atan2( out[k][1], out[k][0] );
				break;
			default:
				for(k=0; k<n; k++) seis[k] = out[k][0];
		}
      //for(k=0; k<n; k+=1000) if(seis[k] != 0.) printf("%d %f\n", k, seis[k]);
   }
   // backward FFT: seis ==> out
   void FFTW_B(int type, float *seis, int n, fftw_complex **in, fftw_complex **out, int *nso, fftw_plan *planF = nullptr) {
      int ns = (int)(log((double)n)/log(2.))+1;
      //if(ns<13) ns = 13;
      ns = (int)pow(2,ns); *nso = ns;
      *in = (fftw_complex *) fftw_malloc ( ns * sizeof(fftw_complex) );//fftw_alloc_complex(ns);
      *out = (fftw_complex *) fftw_malloc ( ns * sizeof(fftw_complex) );//fftw_alloc_complex(ns);
		if( *in==nullptr || *out==nullptr )
			throw ErrorSR::MemError( FuncName, "fftw_malloc failed!");

      //measure plan using the allocated in/out blocks
      //pthread_mutex_lock(&fftlock);
		fftw_plan plan;
		#pragma omp critical(fftw)
      {
      plan = fftw_plan_dft_1d (ns, *in, *out, FFTW_BACKWARD, type); //FFTW_ESTIMATE / FFTW_MEASURE
      if( planF ) *planF = fftw_plan_dft_1d (ns, *out, *in, FFTW_FORWARD, type);
		}
		/* comparing plan to null. Is it a defined behavior
      if( plan==nullptr || (planF && *planF==nullptr) ) {
         fprintf(stderr,"Error(FFTW_B): fftw_plan creation failed!!\n");
         //pthread_mutex_unlock(&fftlock); exit(0);
      }
		*/
      //pthread_mutex_unlock(&fftlock);
      //initialize input array and excute
      memset(*in, 0, ns*sizeof(fftw_complex));
      int k;
      for(k=1; k<n; k++) (*in)[k][0] = seis[k];
      fftw_execute(plan);
		if( !planF ) fftw_free(*in);
      //cleanup
      //pthread_mutex_lock(&fftlock);
		#pragma omp critical(fftw)
      {
      fftw_destroy_plan(plan);
		}
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
   // lowpass taper (cosine)
   void cosTaperR( double f3, double f4, double dom, int nk, fftw_complex *sf, int norderT ) {
      int i = (int)ceil(f3/dom);
		double fwidth = f4-f3, ftmp = M_PI / fwidth;
      for(double f=i*dom; f<f4; i++, f+=dom) {
         double s = ( 1. + cos(ftmp*(f3-f)) ) * 0.5, ss = s;
         for(int io=1; io<norderT; io++) ss *= s;
			sf[i][0] *= ss;
			sf[i][1] *= ss;
      }
      for(;i<nk;i++) {
         sf[i][0] = 0.;
         sf[i][1] = 0.;
      }

      return;
   }
   // lowpass taper (butterworth)
   void btwTaperR( double fc, int norder, double dom, int nk, fftw_complex *sf, int norderT ) {
		// compute stopband freqs
		const double toler = 1.0e-3;
		double fl = fc * pow(1./((1-toler)*(1-toler)) - 1., 0.5/norder);
		double fh = fc * pow(1./(toler*toler) - 1., 0.5/norder);
		// apply taper
		int index = 2 * norder, i = (int)ceil(fl/dom);
      for(double f=i*dom; f<fh&&i<nk; i++, f+=dom) {
			double s = 1. / sqrt(1. + pow(f/fc, index)), ss = s;
         for(int io=1; io<norderT; io++) ss *= s;
         sf[i][0] *= ss;
         sf[i][1] *= ss;
      }
      for(;i<nk;i++) {
         sf[i][0] = 0.;
         sf[i][1] = 0.;
      }

      return;
   }
   // highpass taper (cosine)
   void cosTaperL( double f1, double f2, double dom, int nk, fftw_complex *sf, int norderT ) {
		int i; double f;
		for( i=0, f=0.; f<f1; i++, f+=dom) {
			sf[i][0] = 0.;
			sf[i][1] = 0.;
		}
		double fwidth = f2-f1, ftmp = M_PI / fwidth;
		for(; f<f2; i++, f+=dom) {
			double s = (1. + cos(ftmp*(f2-f))) * 0.5, ss = s;
         for(int io=1; io<norderT; io++) ss *= s;
         sf[i][0] *= ss;
         sf[i][1] *= ss;
		}
      return;
	}
   // highpass taper (butterworth)
   void btwTaperL( double fc, int norder, double dom, int nk, fftw_complex *sf, int norderT ) {
		// compute stopband freqs
		const double toler = 1.0e-3;
		double fl = fc / pow(1./(toler*toler) - 1., 0.5/norder);
		double fh = fc / pow(1./((1-toler)*(1-toler)) - 1., 0.5/norder);
		// apply taper
		int i; double f;
		for( i=0, f=0.; f<fl; i++, f+=dom) {
			sf[i][0] = 0.;
			sf[i][1] = 0.;
		}
		int index = 2 * norder;
      for(; f<fh&&i<nk; i++, f+=dom) {
			double s = 1. / sqrt(1. + pow(fc/f, index)), ss = s;
         for(int io=1; io<norderT; io++) ss *= s;
         sf[i][0] *= ss;
         sf[i][1] *= ss;
      }
      return;
	}
   // bandpass taper (cosine)
   void cosTaperB( double f1, double f2, double f3, double f4, double dom, int nk, fftw_complex *sf, int norderT ) {
		cosTaperL( f1, f2, dom, nk, sf, norderT );
		cosTaperR( f3, f4, dom, nk, sf, norderT );
      return;
   }
   // bandpass taper (butterworth)
   void btwTaperB( double fcL, double fcR, int norder, double dom, int nk, fftw_complex *sf, int norderT ) {
		//btwTaperL( fcL, nL, dom, nk, sf );	// wrong
		//btwTaperR( fcR, nR, dom, nk, sf );	// and wrong!
		// apply taper
		int i; double f;
		int index = 2 * norder;
		for( i=0, f=0.; i<nk; i++, f+=dom) {
			double s = 1. / sqrt(1. + pow((f*f-fcR*fcL)/((fcR-fcL)*f), index)), ss = s;
         for(int io=1; io<norderT; io++) ss *= s;
         sf[i][0] *= ss;
         sf[i][1] *= ss;
      }
      return;
   }
   // gaussian taper
   void gauTaper( double fcenter, double fhlen, double dom, int nk, fftw_complex *sf, int norderT ) {
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
			float ss = gauamp;
         for(int io=1; io<norderT; io++) ss *= gauamp;
			sf[i][0] *= ss;
			sf[i][1] *= ss;
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

	float FillGap( float *pbeg, float *pend, float mean1, float mean2, float std, int hlen, int step ) {
		// defube random number generator
		unsigned timeseed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator (timeseed);
		std::uniform_real_distribution<float> distribution(-std, std);
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
		// shift generated radom numbers to mean1 - mean2
		slope = (mean2 - mean1) / (pend - pbeg);
		for(p=pbeg; p<pend; p++) {
			*p += mean1 + (p-pbeg) * slope;
		}

   }

/*
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
*/

   /* divide the current spectrum by the input freq-amp and freq-phase array */
   void fDiv(double dom, int nsig, fftw_complex *sf, double *freq, double *amp, double *pha, int ntra) {
      int isig, itra = 1;
      double f, ampcur, phacur, realsig, imagsig;
      double sintmp, costmp;
      for(isig=(int)ceil(freq[0]/dom); isig<nsig; isig++) {
         f = isig*dom;
         while(f > freq[itra]) {
	    itra++;
	    if(itra >= ntra) break;
         }
         //interpolate to get current amp and pha
         sintmp = (f-freq[itra-1]) / (freq[itra]-freq[itra-1]);
         ampcur = amp[itra-1] + (amp[itra]-amp[itra-1]) * sintmp;
         phacur = pha[itra-1] + (pha[itra]-pha[itra-1]) * sintmp;
         //divide sf by (ampcur, phacur)
         realsig = sf[isig][0]; imagsig = sf[isig][1];
         sintmp = sin(phacur); costmp = cos(phacur);
         sf[isig][0] = (realsig*costmp - imagsig*sintmp) / ampcur;
         sf[isig][1] = (imagsig*costmp + realsig*sintmp) / ampcur;
      }
   }

   void FDivide (double f1, double f2, double f3, double f4, double *freq, double *amp, double *pha, int nf, SAC_HD &shd, float *seis_in) {
      double dt = shd.delta;
      int n = shd.npts;
      if(f4 > 0.5/dt)
			throw ErrorSR::BadParam( FuncName, "filter band out of range" );
      fftw_plan plan1;
      //backward FFT: s ==> sf
      int ns;
      fftw_complex *s, *sf;
      FFTW_B(FFTW_ESTIMATE, seis_in, n, &s, &sf, &ns, &plan1);
      //make tapering
      int nk = ns/2+1;
      double dom = 1./dt/ns;
      fDiv(dom, nk, sf, freq, amp, pha, nf);
      cosTaperB( f1, f2, f3, f4, dom, nk, sf, 1 );

      //forward FFT: sf ==> s
      //fftw_plan plan2;
      FFTW_F(plan1, s, seis_in, n);
      fftw_free(s); fftw_free(sf);

      //forming final result
      int k;
      float ftmp = 2./(ns*dt);
		for(k=0; k<n; k++) {
			//if( seis_in[k]==0 ) seis_out[k] = 0.;
			//else 
			seis_in[k] *= ftmp;
		}
   }

   int Jday ( int y, int m, int d ) {
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

	#ifdef DISAZI_H
	void ComputeDisAzi( SAC_HD& shd ) {
		try {
			Path<float> path( shd.evlo, shd.evla, shd.stlo, shd.stla );
			shd.dist = path.Dist(); shd.az = path.Azi1(); shd.baz = path.Azi2();
		} catch (const std::exception& e) {}
	}
	#else
	void ComputeDisAzi( SAC_HD& shd ) { throw ErrorSR::UndefMethod( FuncName, "DisAzi.h not included!" ); }
	#endif
};


/* system tools */
#include <sys/types.h>
#include <fts.h>
#include <fnmatch.h>
//#define BLKSIZE 1024
namespace System {
   /* --------------------- Delete file or directory ---------------------------- */
   //#define _XOPEN_SOURCE 500
   //#include <ftw.h>
   //#include <unistd.h>
   int fdprompt = -10;
   void fRemove (const char *fname) {
      //cerr<<"Removing "<<fname<<endl;
      if( remove( fname ) == 0 ) return; //succeed
      int ersv = errno;
      if( ersv == ENOENT ) return; //file not exists
      if( fdprompt>0 ) return; fdprompt++;
      perror("### Warning: Deleting failed"); //failed. prompt to continue
      //TimedContinue(10);
   }

   int Unlink(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf) {
      fRemove((char *)fpath);
      return 0;
   }
//   int dRemove(const char *dirname) {
//      return nftw(dirname, Unlink, 64, FTW_DEPTH | FTW_PHYS);
//   }


   /* ------------------------ Listing (wildcards matching) ----------------------------- */
   //int namecmp(const FTSENT **f1, const FTSENT **f2) { return strcmp((*f1)->fts_name, (*f2)->fts_name); }
   bool List(const char *dir, const char *pattern, int type, std::vector<std::string> &filelist) {
      /* type value decides how sub-directories are handdled
      0: list files in the root dir only
      1: list files in the root dir with dir paths
      2: list all file names
      3: list all files with dir paths */
      //*nfile = 0;
		if( type>3 || type<0 ) {
			//std::cerr<<"ERROR(List): Unknow list type: "<<type<<std::endl;
			//return false;
			throw ErrorSR::BadParam( FuncName, "Unknow list type = "+type );
		}
		FTS *tree;
		FTSENT *file;
		char *dirlist[] = { (char *)dir, nullptr }; //may send in multiple dirs
		//get handle of the file hierarchy; FTS_LOGICAL follows symbolic links and detects cycles.
		//replace '0' with 'namecmp' to sort files by name
		tree = fts_open(dirlist, FTS_LOGICAL | FTS_NOSTAT | FTS_NOCHDIR, 0);
      if (tree == nullptr) perror("fts_open");

      //char *sblk = nullptr;
      //int sleng = 0, bsize = 0;
      //empty the input filelist
      filelist.clear();
      int outflag = 1; // if listing within current directory
      if( type<2 && strcmp(dir, ".")==0 ) outflag=0; // path will not be printed
      //ignores '.' and '..' as FTS_SEEDOT is not set
		while ((file = fts_read(tree))) {
			std::string strtmp;
			switch (file->fts_info) { //current node
				case FTS_DNR: // is a non-readable dir
				case FTS_ERR: // has common errors
				case FTS_NS: // has no stat info
				case FTS_DC: // causes cycle
					//strtmp = "SacSystem::List: "+std::string(file->fts_path)+" causes cycle!";
					//perror( strtmp.c_str() );
				case FTS_DP: // is a post-order dir
					continue; //skip all above cases

				case FTS_D: // is a directory
					if(file->fts_level>0) switch(type) {
						case 0:
							fts_set(tree, file, FTS_SKIP); //no descend
							continue; // and skip
						case 1:
							fts_set(tree, file, FTS_SKIP); //no descend
							break; // and stop switch
						case 2:
							continue; //skip directories
						case 3:;
					}
			}

			if (fnmatch(pattern, file->fts_name, FNM_PERIOD) == 0) {
				/*
					if( sleng > bsize-PLENMAX ) {
					bsize += BLKSIZE;
					sblk = (char *) realloc (sblk, bsize * sizeof(char));
					}
	    if(outflag) sleng += sprintf(&sblk[sleng], "%s\n", file->fts_path);
	    else sleng += sprintf(&sblk[sleng], "%s\n", file->fts_name);
	    */
	    if(outflag) filelist.push_back(file->fts_path);
	    else filelist.push_back(file->fts_name);
	    //*nfile = *nfile+1;
	 }
      }

      if (errno != 0) perror("fts_read");
      if (fts_close(tree) < 0) perror("fts_close");
      if( filelist.size() == 0 ) return false;
      return true;
      //return sblk;
   }

};


//extern MyLogger logger;
/* ---------------------------------------- constructors and operators ---------------------------------------- */
/* default constructor */
SacRec::SacRec( std::ostream& reportin )
 : sig(nullptr), shd(sac_null),
	report(&reportin),
	pimpl(new SRimpl() ) {
}

/* constructor with sac file name */
SacRec::SacRec( const std::string& fnamein, std::ostream& reportin )
 : sig(nullptr), shd(sac_null),
	report(&reportin), fname(fnamein),
	pimpl(new SRimpl() ) {
}

/* constructor with initial signal npts */
SacRec::SacRec( const size_t npts, std::ostream& reportin )
	: SacRec( reportin ) {
	ResizeSig(npts);
}

/* copy constructor */
SacRec::SacRec( const SacRec& recin )
 : fname(recin.fname), report(recin.report),
	shd(recin.shd), sig(new float[recin.shd.npts]), pimpl( new SRimpl(*(recin.pimpl)) ) { 
	if( ! sig )
		throw ErrorSR::MemError( FuncName, "new failed!");
   std::copy(recin.sig.get(), recin.sig.get()+recin.shd.npts, sig.get()); 
}

/* move constructor */
SacRec::SacRec( SacRec&& recin )
 : fname(std::move(recin.fname)), report(recin.report),
	shd(recin.shd), pimpl( std::move(recin.pimpl) ) {
   if ( recin.sig ) {
      sig = std::move(recin.sig);
		//recin.sig.reset(nullptr);
	}
	recin.shd = sac_null;
	recin.report = &(std::cerr);
}

/* assignment operator */
SacRec& SacRec::operator= ( const SacRec& recin ) { 
   pimpl.reset( new SRimpl(*(recin.pimpl)) );
   fname = recin.fname; report = recin.report;
	shd = recin.shd;
   int npts=recin.shd.npts; 
	if( npts > 0 ) {
		sig.reset(new float[npts]);
		if( ! sig )
			throw ErrorSR::MemError( FuncName, "new failed!");
		std::copy(recin.sig.get(), recin.sig.get()+npts, sig.get());
	} else if(sig) {
		sig.reset();
	}
	return *this;
}

/* move ass operator */
SacRec& SacRec::operator= ( SacRec&& recin ) {
	//logger.Hold( WARNING, "In move ass ", FuncName );
	if( this == &recin ) return *this;
   pimpl = std::move(recin.pimpl);
   fname = std::move(recin.fname);
	report = recin.report;	recin.report = &(std::cerr);
   shd = recin.shd; recin.shd = sac_null;
	if ( recin.sig ) {
		sig = std::move(recin.sig);
		recin.sig.reset(nullptr);
	}
	//recin.~SacRec();
	return *this;
}

/* destructor */
SacRec::~SacRec() {}

/* assignment operator */
void SacRec::MutateAs ( const SacRec& recin ) { 
   pimpl.reset( new SRimpl(*(recin.pimpl)) );
   //fname = recin.fname; 
	report = recin.report;
	shd = recin.shd;
	fname = recin.fname;
   sig.reset(new float[recin.shd.npts]);
	if( ! sig )
		throw ErrorSR::MemError( FuncName, "new failed!");
}


/* ------------------------------- memory consumed ------------------------------- */
float SacRec::MemConsumed() const {
	int npts = std::max(0, shd.npts);
	return ( sizeof(SRimpl)+sizeof(SacRec)+fname.size()+npts*sizeof(float) ) / 1048576.;
}

/* ---------------------------------------- sac IO ---------------------------------------- */
/* load sac header from file 'fname' */
void SacRec::LoadHD () {
   std::ifstream fsac(fname.c_str());
   if( ! fsac )
		throw ErrorSR::BadFile( FuncName, "reading from " + fname );
   //if( SHDMap.empty() ) pimpl->CreateSHDMap();
   //pthread_mutex_lock(&fiolock);
	size_t rdsize = sizeof(SAC_HD);
   fsac.read( reinterpret_cast<char *>(&shd), rdsize );
	if( fsac.gcount() != rdsize )
		throw ErrorSR::BadFile( FuncName, "failed to retrieve sac header from " + fname );
   fsac.close();
   //pthread_mutex_unlock(&fiolock);
}

/* load sac header+signal from file 'fname', memory is allocated on heap */
void SacRec::Load () {
   //if( SHDMap.empty() ) pimpl->CreateSHDMap();
   //sig = std::make_shared<float>( new float[shd.npts*sizeof(float)] );
	// check input file
   std::ifstream fsac(fname.c_str());
   if( ! fsac )
		throw ErrorSR::BadFile( FuncName, "reading from " + fname );
	// read from fin
	size_t rdsize = sizeof(SAC_HD);
	#pragma omp critical(sacIO)
	{
   fsac.read( reinterpret_cast<char *>(&shd), rdsize );
	if( fsac.gcount() != rdsize )
		throw ErrorSR::BadFile( FuncName, "failed to retrieve sac header from " + fname );
	// allocate memory for sac signal
   sig.reset(new float[shd.npts]);
	if( ! sig )
		throw ErrorSR::MemError( FuncName, "new failed!");
   float* sigsac = sig.get();
	rdsize = sizeof(float) * shd.npts;
   fsac.read( reinterpret_cast<char *>(sigsac), rdsize );
	if( fsac.gcount() != rdsize )
		throw ErrorSR::BadFile( FuncName, "failed to extract sac sig from " + fname );
   fsac.close();
	} // critical

   /* calculate t0 */
   char koo[9];
   for ( int i=0; i<8; i++ ) koo[i] = shd.ko[i]; koo[8] = 0;
   float fes;
   int eh, em;
   sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);
   shd.o = shd.b + (shd.nzhour-eh)*3600. + (shd.nzmin-em)*60. + shd.nzsec-fes + shd.nzmsec*.001;
}

/* write to file '*outfname' */
void SacRec::WriteHD (const std::string& outfname) {
   // open file with ios::in|ios::out (no ios::trunc) to avoid deleting the contents
   //std::fstream fsac(outfname, std::ios::in | std::ios::out);
   std::fstream fsac(outfname);	// fstream defaults as std::ios::in|stdLLios::out
   if( ! fsac )
		throw ErrorSR::BadFile( FuncName, "writing to " + std::string(outfname) );
   /* update header */
   shd.iftype = (int)ITIME;
   shd.leven = (int)TRUE;
   shd.lovrok = (int)TRUE;
   shd.internal4 = 6L;

   /* check and re-format header time if necessary */
   UpdateTime();

   //pthread_mutex_lock(&fiolock);
   fsac.write( reinterpret_cast<char *>(&shd), sizeof(SAC_HD) );
   //pthread_mutex_unlock(&fiolock);

   fsac.close();
}

void SacRec::updateDeps() {
   /* search min and max amplitude */
	float *sigsac = sig.get();
   shd.depmin = sigsac[0];
   shd.depmax = sigsac[0];
   for ( int i = 0; i < shd.npts ; i++ ) {
       if ( shd.depmin > sigsac[i] ) shd.depmin = sigsac[i];
       else if ( shd.depmax < sigsac[i] ) shd.depmax = sigsac[i];
   }
}

void SacRec::Write (const std::string& outfname) {
   /* check if signal is loaded */
   if( ! sig )
		throw ErrorSR::EmptySig( FuncName, "writing to " + outfname );
   /* open file */
   std::ofstream fsac(outfname);
   if( ! fsac )
		throw ErrorSR::BadFile( FuncName, "writing to " + outfname );
   /* update header */
   shd.iftype = (int)ITIME;
   shd.leven = (int)TRUE;
   shd.lovrok = (int)TRUE;
   shd.internal4 = 6L;

	// update depmin and depmax
	updateDeps();
   /* check and re-format header time if necessary */
   UpdateTime();

   //pthread_mutex_lock(&fiolock);
	#pragma omp critical(sacIO)
	{
   fsac.write( reinterpret_cast<char *>(&shd), sizeof(SAC_HD) );
   fsac.write( reinterpret_cast<char *>(sig.get()), sizeof(float)*shd.npts );
   fsac.close();
	}
}

void SacRec::LoadTXT( const std::string& fname ) {
	
   std::ifstream fin(fname);
   if( ! fin )
		throw ErrorSR::BadFile( FuncName, "reading from " + fname );
	// read from fin
	std::vector< std::array<float, 2> > dataV;
	float delta = -123.;
	#pragma omp critical(sacIO)
	{
	float xold = -123.;
	for(std::string line; std::getline(fin, line); ) {
		std::stringstream ss(line);
		float x, y; 
		if( !(ss >> x >> y) ) continue;
		dataV.push_back( {x, y} );
		if( delta!=-123. ) {
			float factor = (x - xold) / delta;
			if( factor != (int)factor )
				throw ErrorSR::HeaderMismatch( FuncName, "irregular x" );
		} else if( xold!=-123. ) { delta = x - xold; }
		xold = x;
	}
	}
	// write header
	shd.delta = delta;
	shd.b = dataV.front()[0]; shd.e = dataV.back()[0];
	shd.npts = (shd.e - shd.b) / delta + 1;
	// allocate memory for sac signal
   sig.reset(new float[shd.npts]());
	if( ! sig )
		throw ErrorSR::MemError( FuncName, "new failed!");
   float* sigsac = sig.get();
	// fill signal
	for( const auto& dp : dataV ) {
		int ix = (dp[0] - shd.b) / delta;
		sigsac[ix] = dp[1];
	}
}

void SacRec::Dump( const std::string fname ) const {
	if( !sig || shd.npts<=0 )
		throw ErrorSR::EmptySig(FuncName);
	bool tofile = !fname.empty();
	std::ofstream fout(fname);
	if( tofile && !fout )
		throw ErrorSR::BadFile( FuncName, "writing to " + fname );
	std::ostream& sout = tofile ? fout : std::cout;
	float b = shd.b, dt = shd.delta, *sigsac = sig.get();
	for( int i=0; i<shd.npts; i++ )
		sout<<b+i*dt<<" "<<sigsac[i]<<"\n";
}

/* ------------------------------- header operations ------------------------------- */
void SacRec::PrintHD( const std::string field, std::ostream &o ) const {
	shd.StreamTo(field, o); 
}
void SacRec::DumpHD( const std::string fname ) const {
	bool tofile = !fname.empty();
	std::ofstream fout(fname);
	if( tofile && !fout )
		throw ErrorSR::BadFile( FuncName, "writing to " + fname );
	std::ostream& sout = tofile ? fout : std::cout;

	sout << shd;
}

/*
int read_rec(int rec_flag, char *fname, int len, int *rec_b, int *rec_e, int *nrec) {
   FILE *frec;
   int irec;
   if( rec_flag ) {
      if((frec = fopen(fname,"r")) == nullptr) return 0;
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
void SacRec::Mul( const float mul ) {
   if( !sig )
		throw ErrorSR::EmptySig(FuncName);
	float *sigsac = sig.get();
   for(int i=0; i<shd.npts; i++) sigsac[i] *= mul;
}

void SacRec::Addf( const SacRec& sac2 ) {
   if( !sig || !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );
	//const auto& sig2 = sac2.sig;
	float *sigsac = sig.get(), *sigsac2 = sac2.sig.get();
   for(int i=0; i<shd.npts; i++) 
		sigsac[i] += sigsac2[i];
}

void SacRec::Subf( const SacRec& sac2 ) {
   if( !sig || !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );
	//const auto& sig2 = sac2.sig;
	float *sigsac = sig.get(), *sigsac2 = sac2.sig.get();
   for(int i=0; i<shd.npts; i++) 
		sigsac[i] -= sigsac2[i];
}

void SacRec::Divf( const SacRec& sac2 ) {
   if( !sig || !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );
	//const auto& sig2 = sac2.sig;
	float *sigsac = sig.get(), *sigsac2 = sac2.sig.get();
   for(int i=0; i<shd.npts; i++) 
		if( sigsac2[i]!=0. ) sigsac[i] /= sigsac2[i];
}

float SacRec::SNR( const float tsignall, const float tsignalh, const float tnoisel, const float tnoiseh ) const {
	float tmin, tmax, min, max;
	MinMax( tsignall, tsignalh, tmin, min, tmax, max );
	float Asig = std::max(max, -min);
	float RMSnoise = RMSAvg(tnoisel, tnoiseh);
	float SNR = Asig / RMSnoise;
	//shd.user1 = SNR;
	return SNR;
}

void SacRec::PullUpTo( const SacRec& sac2 ) {
   if( !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( !sig )
		*this = sac2;
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );
	//const auto& sig2 = sac2.sig;
	float *sigsac = sig.get(), *sigsac2 = sac2.sig.get();
   for(int i=0; i<shd.npts; i++) 
		if( sigsac2[i] > sigsac[i] ) sigsac[i] = sigsac2[i];
		//sigsac[i] = std::max(sigsac[i], sigsac2[i]);
}

void SacRec::Reverse( SacRec& sac2 ) {
   if( !sig )
		throw ErrorSR::EmptySig(FuncName);

	if( &sac2 == this )
		throw ErrorSR::BadParam(FuncName, "reverse in place is not allowed!");

	// initialize sac2
   sac2.MutateAs(*this);

	float *sigsac = sig.get(), *sigsac2 = sac2.sig.get();
   for(int i=0; i<shd.npts; i++) 
		sigsac2[shd.npts-i-1] = sigsac[i];
}



float SacRec::Dis() const {
	auto shdl = shd;
	if( shdl.dist == NaN ) pimpl->ComputeDisAzi( shdl );
	return shdl.dist;
}
float SacRec::Dis() {
	if( shd.dist == NaN ) pimpl->ComputeDisAzi( shd );
	return shd.dist;
}
float SacRec::Azi() const {
	auto shdl = shd;
	if( shdl.az == NaN ) pimpl->ComputeDisAzi( shdl );
	return shdl.az;
}
float SacRec::Azi() {
	if( shd.az == NaN ) pimpl->ComputeDisAzi( shd );
	return shd.az;
}

double SacRec::AbsTime () {
   //if( ! sig ) return -1.;
   //if( shd == sac_null ) return -1.; // operator== not defined yet
   if( shd.npts <= 0 ) return -1.;
   if( shd.nzjday == NaN || shd.nzyear == NaN || shd.nzhour == NaN ||
       shd.nzmin == NaN || shd.nzsec == NaN || shd.nzmsec == NaN ) return -1;
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
/*
void SacRec::MinMax (int& imin, int& imax) const {
	float min, max;
	float *sigsac = sig.get();
	min = max = sigsac[0];
   imin = imax = 0;
	for( int i=0; i<shd.npts; i++ ) {
       if ( min > sigsac[i] ) { min = sigsac[i]; imin = i; }
       else if ( max < sigsac[i] ) { max = sigsac[i]; imax = i; }
   }
}
*/
void SacRec::MinMax (int& imin, int& imax, float tbegin, float tend) const {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);
	// correct time window if necessary
	if( tbegin==NaN || tbegin<shd.b ) tbegin = shd.b;
	if( tend==NaN || tend>shd.e ) tend = shd.e;
	// indexing range
	const int il = Index(tbegin), ih = Index(tend)+1;
   imin = imax = il;
	// initialize min, max
	float *sigsac = sig.get();
	float min, max; min = max = sigsac[il];
   for ( int i = il+1; i < ih; i++ ) {
       if ( min > sigsac[i] ) { min = sigsac[i]; imin = i; }
       else if ( max < sigsac[i] ) { max = sigsac[i]; imax = i; }
   }
}

void SacRec::MinMax (float tbegin, float tend, float& tmin, float& min, float& tmax, float& max) const {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);
	// correct time window if necessary
	if( tbegin == NaN ) tbegin = shd.b;
	if( tend == NaN ) tend = shd.e;
	// indexing range
	const int il = Index(tbegin), ih = Index(tend)+1;
   int imin = il, imax = il;
	float *sigsac = sig.get();
	min = max = sigsac[il];
   for ( int i = il+1; i < ih ; i++ ) {
       if ( min > sigsac[i] ) { min = sigsac[i]; imin = i; }
       else if ( max < sigsac[i] ) { max = sigsac[i]; imax = i; }
   }
	// time results
   tmin = shd.b + imin*shd.delta;
   tmax = shd.b + imax*shd.delta;
}

/* compute the root-mean-square average in a given window */
float SacRec::RMSAvg ( float tbegin, float tend, int step ) const {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);
	/* decide avg window */
   float maxfloat = std::numeric_limits<float>::max();
   int neff = 0, ibeg = Index(tbegin), iend = Index(tend)+1;
	/*
	if( ibeg < 0 ) {
		ibeg = 0;
		(*report) << "RMSAvg: tbegin out of range!" << std::endl;
	}
	if( iend > shd.npts ) {
		iend = shd.npts;
		(*report) << "RMSAvg: tend out of range!" << std::endl;
	}
	if( ibeg >= iend )
		throw ErrorSR::BadParam(FuncName, "0 window length");
	*/
	if( ibeg<0 || iend>shd.npts || ibeg >= iend )
		throw ErrorSR::BadParam(FuncName, "invalid time window input");

   float rms = 0.;
	float *sigsac = sig.get();
   for ( int i = ibeg; i < iend ; i+=step ) {
      if( sigsac[i] >= maxfloat ) continue;
      rms += sigsac[i] * sigsac[i];
      neff++;
   }
   rms = sqrt(rms/(neff-1));
	return rms;
}


/* compute accurate time of the peak (fit with a parabola) */
void SacRec::Peak(float& tpeak, float& apeak, const float tbegin, const float tend) const {
	// find isig closest to the peak
   int imin, imax;
   MinMax(imin, imax, tbegin, tend );
   float* sigsac = sig.get();
	if( fabs(sigsac[imin])>fabs(sigsac[imax]) ) imax = imin;
	// fit a parabola to get the accurate peak time
#ifndef PARABOLA_H
	tpeak = Time(imax);
	apeak = sigsac[imax];
#else
   if( imax==0 || imax==shd.npts-1 ) {
      tpeak = X(imax);
		apeak = sigsac[imax];
   } else {
      PointC p1(X(imax-1), sigsac[imax-1]);
      PointC p2(X(imax),   sigsac[imax]  );
      PointC p3(X(imax+1), sigsac[imax+1]);
      auto v = Parabola( p1, p2, p3 ).Vertex();
		tpeak = v.x; apeak = v.y;
   }
#endif
}

/* compute accurate sig value at a given time (fit with a parabola) */
float SacRec::Sig( float time ) const {
	int itime = Index(time);
#ifndef PARABOLA_H
	return sig.at(itime);
#else
	if( itime == 0 ) itime += 1;
	else if( itime == shd.npts-1 ) itime -= 1;
	float* sigsac = sig.get();
	PointC p1(X(itime-1), sigsac[itime-1]);
	PointC p2(X(itime),   sigsac[itime]  );
	PointC p3(X(itime+1), sigsac[itime+1]);
	//std::cerr<<p1<<"\n"<<p2<<"\n"<<p3<<"\nresult: "<<time<<" "<<Parabola( p1, p2, p3 )(time)<<std::endl;
	return Parabola( p1, p2, p3 )(time);
#endif
}


/* smoothing (running average) */
void SacRec::Smooth( float timehlen, SacRec& sacout ) const {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);

	/* resize sacout.sig */
	sacout.shd = shd;
	sacout.sig.reset( new float[shd.npts] );
	if( ! sig )
		throw ErrorSR::MemError( FuncName, "new failed!");

	/* copy and compute abs of the signal into sigw */
   int i, j, wb, we, n = shd.npts;
   float wsum, dt = shd.delta;
   int half_l = (int)floor(timehlen/dt+0.5);
   float *sigsac = sig.get(), *sigout = sacout.sig.get();
	//float sigw[n];
	std::unique_ptr<float> sigw_p(new float[n]);
   float *sigw = sigw_p.get(); 
	for( i=0; i<n; i++ ) sigw[i] = fabs(sigsac[i]);

	/* */
   if(half_l*2>n-1) half_l = (n-1)/2;
   for(i=0,wsum=0.;i<=half_l;i++) wsum += sigw[i];
   wb = 0; we = i;
   for(i=1;i<=half_l;i++,we++) {
		sigout[i-1] = (double)wsum/we;
      wsum += sigw[we];
   }
   for(j=we;i<n-half_l;i++,wb++,we++) {
		//if( i>80000/dt && i<82000/dt ) std::cerr<<(i-1)*dt<<" "<<sig[i-1]<<" "<<wsum<<" / "<<j<<std::endl;
		sigout[i-1] = (double)wsum/j;
      wsum += ( sigw[we] - sigw[wb] );
   }
   for(;i<n;i++,wb++) {
		sigout[i-1] = (double)wsum/(we-wb);
      wsum -= sigw[wb];
   }
   if(wsum>1.e-15)
		sigout[n-1] = (double)wsum/(we-wb);

}


/* compute mean and std in a given window */
bool SacRec::Mean ( float tbegin, float tend, int step, float& mean ) const {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);

	// data window
   float maxfloat = std::numeric_limits<float>::max();
   int ibeg = Index(tbegin), iend = Index(tend)+1;
	if( ibeg<0 || iend>shd.npts ) return false;

	//compute mean on valid data points
	float *sigsac = sig.get();
	mean = 0.; int neff = 0;
   for ( int i = ibeg; i < iend ; i+=step ) {
      if( sigsac[i] >= maxfloat ) continue;
		mean += sigsac[i];
		neff++;
   }
	if( neff == 0 ) return false;
	mean /= neff;
}

bool SacRec::MeanStd ( float tbegin, float tend, int step, float& mean, float& std ) const {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);

	// store valid data points
   float maxfloat = std::numeric_limits<float>::max();
   int ibeg = Index(tbegin), iend = Index(tend)+1;
	if( ibeg<0 || iend>shd.npts ) return false;

	std::vector<float> dataV;
	float *sigsac = sig.get();
   for ( int i = ibeg; i < iend ; i+=step ) {
      if( sigsac[i] >= maxfloat ) continue;
      dataV.push_back(sigsac[i]);
   }
	if( dataV.size() == 0 ) return false;
	//	throw ErrorSR::InsufData(FuncName, "0 valid data point within the given window");

	// mean
   mean = 0.;
	for( int i=0; i<dataV.size(); i++ ) mean += dataV[i];
	mean /= dataV.size();
	// std
	std = 0.;
	for( int i=0; i<dataV.size(); i++ ) {
		float ftmp = dataV[i] - mean;
		std += ftmp * ftmp;
	}
   std = sqrt(std/(dataV.size()-1));

	return true;
}


/* phase wrap and unwrap */
void SacRec::Wrap() {
	Transform( [](float& val){ val -= nint(val/twopi) *twopi;	} );
}
void SacRec::Unwrap() {
	if( ! sig )
		throw ErrorSR::EmptySig(FuncName);
	float *sigph = sig.get();
	for(int i=1; i<shd.npts; i++)
		sigph[i] -= nint( (sigph[i]-sigph[i-1])/twopi ) * twopi;
}

/* performs integration in time domain using the trapezoidal rule */
void SacRec::IntegrateT( SacRec& sac_out ) const {
	if( &sac_out != this ) sac_out = *this;
	float dt = shd.delta, hdt = dt * 0.5;
	float sum = -hdt * sig[0];
	sac_out.Transform( [&](float& val) {
		float fadd = val * hdt;
		sum += fadd;
		val = sum;
		sum += fadd;
	} );
}

/* performs integration in the frequency domain (omega arithmetic) */
void SacRec::Integrate( SacRec& sac_out ) const {
	// FFT
	SacRec sac_am, sac_ph;
	ToAmPh( sac_am, sac_ph );
	// shift phase by -pi/2 (FFT is backward!)
	float HLF_PI = M_PI * 0.5, TWO_PI = M_PI * 2.;
	sac_ph.Transform( [&](float& val) {
		val += HLF_PI;
		if( val >= M_PI ) val -= TWO_PI;
	} );
	// lowpass filtering at 1000 sec
	float fl = 0.001; size_t ifl = sac_am.Index(fl);
	sac_am.Transform( [](float& val) { 
		val = 0.; 
	}, 1, ifl );
	// and divide amplitude by omega
	float domega = sac_am.shd.delta * 2. * M_PI, omega = ifl * domega;
	sac_am.Transform( [&](float& val) {
		val /= omega;
		omega += domega;
	}, fl );
	// IFFT
	sac_out.shd = shd;
	sac_out.FromAmPh( sac_am, sac_ph );
}

/* performs differentiation in the frequency domain (omega arithmetic) */
void SacRec::Differentiate( SacRec& sac_out ) const {
	// FFT
	SacRec sac_am, sac_ph;
	ToAmPh( sac_am, sac_ph );
	// shift phase by pi/2 (FFT is backward!)
	float HLF_PI = M_PI * 0.5, TWO_PI = M_PI * 2.;
	sac_ph.Transform( [&](float& val) {
		val -= HLF_PI;
		if( val < -M_PI ) val += TWO_PI;
	} );
	// multiply amplitude by omega
	float domega = sac_am.shd.delta * 2. * M_PI, omega = 0.;
	sac_am.Transform( [&](float& val) {
		val *= omega;
		omega += domega;
	} );
	// IFFT
	sac_out.shd = shd;
	sac_out.FromAmPh( sac_am, sac_ph );
}


/* ---------------------------------------- time - frequency ---------------------------------------- */
/* convert to amplitude */
void SacRec::ToAmPh( SacRec& sac_am, SacRec& sac_ph, const int nfout ) const {
	if( shd.npts > maxnpts4parallel ) {
		#pragma omp critical(largesig)
		{
		ToAmPh_p(sac_am, sac_ph, nfout);
		}
	} else {
		ToAmPh_p(sac_am, sac_ph, nfout);
	}
}
void SacRec::ToAmPh_p( SacRec& sac_am, SacRec& sac_ph, const int nfout ) const {
   if( ! sig ) 	// check signal
		throw ErrorSR::EmptySig(FuncName);
   if( &sac_am != this ) {	// initialize sac_am if not filtering in place 
      //sac_am = *this;
      sac_am.fname = fname; 
      sac_am.shd = shd;
      sac_am.pimpl.reset( new SRimpl(*(pimpl)) );
   }

   //fftw_plan planF = nullptr;
   //backward FFT: s ==> sf
   int ns, n=shd.npts;
   fftw_complex *s, *sf;
   pimpl->FFTW_B(FFTW_ESTIMATE, &(sig[0]), n, &s, &sf, &ns);

   //forming amplitude spectrum
   int nk = std::max(nfout, ns/2 + 1);
   float *amp = new float[nk], *pha = new float[nk];
	float delta = shd.delta;
	if( amp==nullptr || pha==nullptr )
		throw ErrorSR::MemError( FuncName, "new failed!");
   for(int i=0; i<nk; i++) {
      fftw_complex& cur = sf[i];
      amp[i] = sqrt(cur[0]*cur[0] + cur[1]*cur[1]) * delta;
		pha[i] = atan2(cur[1], cur[0]);
   }
   //fftw_free(s); 
	fftw_free(sf);
   sac_am.sig.reset(amp);
	sac_ph.sig.reset(pha);

   sac_am.shd.npts = nk;
   sac_am.shd.delta = 1./(delta*ns);
   sac_am.shd.b = 0.;
	sac_am.shd.e = sac_am.shd.delta * nk;
	sac_ph.shd = sac_am.shd;
}

/* convert to re & im */
void SacRec::FFT( SacRec& sac_re, SacRec& sac_im, const int nfout ) const {
	if( shd.npts > maxnpts4parallel ) {
		#pragma omp critical(largesig)
		{
		FFT_p(sac_re, sac_im, nfout);
		}
	} else {
		FFT_p(sac_re, sac_im, nfout);
	}
}
void SacRec::FFT_p( SacRec& sac_re, SacRec& sac_im, const int nfout ) const {
   if( ! sig ) 	// check signal
		throw ErrorSR::EmptySig(FuncName);
   if( &sac_re != this ) {	// initialize sac_re if not filtering in place 
      //sac_re = *this;
      sac_re.fname = fname; 
      sac_re.shd = shd;
      sac_re.pimpl.reset( new SRimpl(*(pimpl)) );
   }

   //fftw_plan planF = nullptr;
   //backward FFT: s ==> sf
   int ns, n=shd.npts;
   fftw_complex *s, *sf;
   pimpl->FFTW_B(FFTW_ESTIMATE, &(sig[0]), n, &s, &sf, &ns);

   //forming spectrums
   int nk = std::max(nfout, ns/2 + 1);
   float *re = new float[nk], *im = new float[nk];
	float delta = shd.delta;
	if( re==nullptr || im==nullptr )
		throw ErrorSR::MemError( FuncName, "new failed!");
   for(int i=0; i<nk&&i<ns; i++) {
      fftw_complex& cur = sf[i];
      re[i] = cur[0] * delta;
		im[i] = cur[1] * delta;
   }
   //fftw_free(s); 
	fftw_free(sf);
   sac_re.sig.reset(re);
	sac_im.sig.reset(im);

   sac_re.shd.npts = nk;
   sac_re.shd.delta = 1./(delta*ns);
   sac_re.shd.b = 0.;
	sac_re.shd.e = sac_re.shd.delta * nk;
	sac_im.shd = sac_re.shd;
}

void SacRec::FromAmPh( SacRec& sac_am, SacRec& sac_ph, const short outtype ) {
	if( shd.npts > maxnpts4parallel ) {
		#pragma omp critical(largesig)
		{
		FromAmPh_p(sac_am, sac_ph, outtype);
		}
	} else {
		FromAmPh_p(sac_am, sac_ph, outtype);
	}
}
void SacRec::FromAmPh_p( SacRec& sac_am, SacRec& sac_ph, const short outtype ) {
	/* outtype =
		0: original
		1: hilbert
		2: envelope
		3: signal phase (in time domain)
	*/
	/*
	// check header of current sac
	if( shd.npts==NaN || shd.delta==NaN || shd.b==NaN )
		throw ErrorSR::BadParam(FuncName, "empty header");
	*/
	// check signals
   if( !sac_am.sig || !sac_ph.sig )
      throw ErrorSR::EmptySig(FuncName);
	if( sac_am.shd.delta != sac_ph.shd.delta )
		throw ErrorSR::HeaderMismatch(FuncName, "dt");
	int nk = sac_am.shd.npts;
	if( nk != sac_ph.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(nk)+" - "+std::to_string(sac_ph.shd.npts) );
	// allocate memory
	int ns = (nk-1) * 2;
	fftw_complex *in = (fftw_complex *) fftw_malloc ( ns * sizeof(fftw_complex) );//fftw_alloc_complex(ns);
	fftw_complex *out = (fftw_complex *) fftw_malloc ( ns * sizeof(fftw_complex) );//fftw_alloc_complex(ns);
	if( *in==nullptr || *out==nullptr )
		throw ErrorSR::MemError( FuncName, "fftw_malloc failed!");
	// assign spectrum to the 'in' array
   memset(in, 0, ns*sizeof(fftw_complex));
	float *sigam = sac_am.sig.get(), *sigph = sac_ph.sig.get();
   for(int i=0; i<nk; i++) {
		in[i][0] = sigam[i] * cos(sigph[i]);
		in[i][1] = sigam[i] * sin(sigph[i]);
   }
	// create plan
	fftw_plan plan;
	#pragma omp critical(fftw)
	{
	plan = fftw_plan_dft_1d (ns, in, out, FFTW_FORWARD, FFTW_ESTIMATE); //FFTW_ESTIMATE / FFTW_MEASURE
	}
	// set header
	if( shd.npts==NaN || shd.npts>ns ) shd.npts = ns;	//shd.npts = ns;
	if( shd.delta == NaN ) shd.delta = 1.;
	// run FFTW_F
	sig.reset( new float[shd.npts] );
	if( ! sig )
		throw ErrorSR::MemError( FuncName, "new failed!");
	pimpl->FFTW_F(plan, out, sig.get(), shd.npts, outtype);
	// free memory
	fftw_free(in); fftw_free(out);
	// normalize
	if( outtype != 3 ) {
		float *sacsig = sig.get(), ftmp = 2./(shd.delta*ns);
		for(int i=0; i<shd.npts; i++)	sacsig[i] *= ftmp;
	}
}

/* method that performs different types of filters:
 * type = 0: Lowpass cosine -f3~f4_
 * type = 1: highpass cosine _f1~f2-
 * type = 2: bandpass cosine _f1~f2-f3~f4_
 * type = 3: lowpass butterworth -fc=f3~n=f4_
 * type = 4: highpass butterworth _fc=f1~n=f2-
 * type = 5: bandpass butterworth _fcL=f2~fcR=f3~n=f4_
 * type = 6: gaussian _fc=f2~fhlen=f3_ */
void SacRec::Filter ( double f1, double f2, double f3, double f4, const int type, SacRec& srout, bool zeroPhase ) {
	if( shd.npts > maxnpts4parallel ) {
		#pragma omp critical(largesig)
		{
		Filter_p(f1, f2, f3, f4, type, srout, zeroPhase);
		}
	} else {
		Filter_p(f1, f2, f3, f4, type, srout, zeroPhase);
	}
}
void SacRec::Filter_p ( double f1, double f2, double f3, double f4, const int type, SacRec& srout, bool zeroPhase ) {
   if( ! sig )	// check signal
		throw ErrorSR::EmptySig(FuncName);
   if( &srout != this ) {	// initialize srout if not filtering in place 
      srout = *this;
      /*
      srout.fname = fname; srout.shd = shd; 
      srout.sig.reset( new float[shd.npts] );
      srout.pimpl.reset( new SRimpl(*(pimpl)) );
      */
   }

   double dt = shd.delta;
	// correct f4 for cos filters
   if(type<3 && f4>0.5/dt) {
      //std::cerr<<"Warning(SacRec::Filter): filter band out of range!"<<std::endl;;
      f4 = 0.49999/dt;
   }
   fftw_plan planF;

   // backward FFT: s ==> sf
   int ns, n=shd.npts;

	// run in series when npts is too large to prevent memory problem
   fftw_complex *s, *sf;
   pimpl->FFTW_B(FFTW_ESTIMATE, &(sig[0]), n, &s, &sf, &ns, &planF);

   //make tapering
   int nk = ns/2+1, norder = zeroPhase ? 2 : 1;
   double dom = 1./(dt*ns);
	switch(type) {
		case 0:	// Lowpass cosine
			pimpl->cosTaperR( f3, f4, dom, nk, sf, norder ); break;
		case 1:	// Highpass cosine
			pimpl->cosTaperL( f1, f2, dom, nk, sf, norder ); break;
		case 2:	// Bandpass cosine
			pimpl->cosTaperB( f1, f2, f3, f4, dom, nk, sf, norder ); break;
		case 3:	// Lowpass butterworth
			pimpl->btwTaperR( f3, f4, dom, nk, sf, norder ); break;
		case 4:	// Highpass butterworth
			pimpl->btwTaperL( f1, f2, dom, nk, sf, norder ); break;
		case 5:	// Bandpass butterworth
			pimpl->btwTaperB( f2, f3, f4, dom, nk, sf, norder ); break;
		case 6:	// Gaussian
			pimpl->gauTaper( f2, f3, dom, nk, sf, norder ); break;
		default:
			throw ErrorSR::BadParam( FuncName, "Unknown filter type");
	}
   //if( (f1==-1. || f2==-1.) && (f3>0. && f4>0.) ) pimpl->cosTaperR( f3, f4, dom, nk, sf );
   //else if( f1>=0. && f2>0. && f3>0. && f4>0. ) pimpl->cosTaperB( f1, f2, f3, f4, dom, nk, sf );
   //else if( f1==-1. && f4==-1. ) pimpl->gauTaper( f2, f3, dom, nk, sf );

   //forward FFT: sf ==> s
   pimpl->FFTW_F(planF, s, &(srout.sig[0]), n);
   fftw_free(s); fftw_free(sf);

   //forming final result
   float ftmp = 2./ns;
	float *sigout = srout.sig.get();
   for(int k=0; k<n; k++) sigout[k] *= ftmp;
   //   if( seis_in[k]==0 ) seis_out[k] = 0.;
   //   else seis_out[k] *= ftmp;
}


void SacRec::gauTaper( const float fc, const float fh ) {
	if( !sig || shd.npts<=0 )	// check signal
		throw ErrorSR::EmptySig(FuncName);

	int i;
	float gauamp, dom = shd.delta;
	double f, fstart, fend, fmax = (shd.npts-1)*shd.npts;
	// define effective window for given gaussian halflength
	f = fh * 4.;
	fstart = fc - f;
	fend = fc + f; 
	if( fend > fmax ) fend = fmax;
	float *sigsac = sig.get();
	// cut high frequecies
	if( fstart > shd.b ) {
		for(i=0, f=shd.b; f<fstart; i++, f+=dom)
			sigsac[i] = 0.;
		f = shd.b + i * dom; // correct for round-off error
	}
	else { f = shd.b; i = 0; }
	// apply taper
	float alpha = -0.5/(fh*fh);
	for(; f<fend-1.e-10; i++, f+=dom) {
		gauamp = f - fc;
		gauamp = exp( alpha * gauamp * gauamp );
		sigsac[i] *= gauamp;
	}
	// cut low frequencies
	if( fend < fmax ) {
		f = shd.b + i * dom; // again, correct for round-off
		for(; i<shd.npts; i++) {
			sigsac[i] = 0.;
		}
	}
}

/* cosine tapper */
void SacRec::cosTaperL( const float fl, const float fh ) {
	if( !sig || shd.npts<=0 )	// check signal
		throw ErrorSR::EmptySig(FuncName);

	const float& delta = shd.delta;
	int i;
	float *sigsac = sig.get();
	for( i=0; i<(int)ceil((fl-shd.b)/delta); i++ ) sigsac[i] = 0.;
	float finit = shd.b+i*delta, fwidth = fh-fl, dfh = fh-finit;
	float ftmp = M_PI / fwidth;
   for(; dfh>0; i++, dfh-=delta) {
		float amplif = ( 1. + cos(ftmp*dfh) ) * 0.5;
		sigsac[i] *= amplif;
   }

   return;
}
void SacRec::cosTaperR( const float fl, const float fh ) {
   if( !sig || shd.npts<=0 )	// check signal
		throw ErrorSR::EmptySig(FuncName);

	const float& delta = shd.delta;
   int i = (int)ceil((fl-shd.b)/delta);
	float finit = shd.b+i*delta, fwidth = fh-fl, dfl = finit-fl;
	float ftmp = M_PI / fwidth;
	float *sigsac = sig.get();
   for(; i<shd.npts&&dfl<fwidth; i++, dfl+=delta) {
		float amplif = ( 1. + cos(ftmp*dfl) ) * 0.5;
		sigsac[i] *= amplif;
   }
   for(;i<shd.npts;i++) sigsac[i] = 0.;

   return;
}

void SacRec::Hilbert( SacRec& sacout ) {
   if( !sig || shd.npts<=0 )	// check signal
		throw ErrorSR::EmptySig(FuncName);

	SacRec sac_am, sac_ph;
	this->ToAmPh(sac_am, sac_ph);
	sacout.FromAmPh(sac_am, sac_ph, 1);	// take the imaginary part from IFFT for hilbert
}

void SacRec::Envelope( SacRec& sacout ) {
   if( !sig || shd.npts<=0 )	// check signal
		throw ErrorSR::EmptySig(FuncName);

	SacRec sac_am, sac_ph;
	this->ToAmPh(sac_am, sac_ph);
	sacout.FromAmPh(sac_am, sac_ph, 2);	// take the amplitude from IFFT for envelope
}


/* ---------------------------------------- cut and merge ---------------------------------------- */
void SacRec::cut( float tb, float te, SacRec& sac_result ) {
   int nb = nint( (tb-shd.b) / shd.delta );
   int ne = nint( (te-shd.b) / shd.delta );
   if( nb>=ne || ne<0 || nb>shd.npts )
		throw ErrorSR::BadParam( FuncName, "Invalid nb/ne/npts = " + std::to_string(nb) +
										 "/" + std::to_string(ne) + "/" + std::to_string(shd.npts) );
   int nptsnew = ne - nb + 1;
   float* signew = (float *) calloc ( nptsnew, sizeof(float) );
	if( signew == nullptr )
		throw ErrorSR::MemError( FuncName, "calloc failed!");
   // define start positions
   int inew, iold;
   if( nb < 0 ) { inew = -nb; iold = 0; }
   else { inew = 0; iold = nb; }
   // define copy size
   float nptscpy = std::min( nptsnew - inew, shd.npts - iold - 1 );
   // copy data
   memcpy( &(signew[inew]), &(sig[iold]), nptscpy * sizeof(float) );
   // reset sacT.sig
   sac_result.sig.reset(signew);
   // update shd
   if( &(sac_result) != this ) {
		sac_result.shd = shd;
		sac_result.fname = fname;
	}
   sac_result.shd.b += nb * shd.delta;
	sac_result.shd.e = sac_result.shd.b + (nptsnew-1) * shd.delta;
   sac_result.shd.npts = nptsnew;
}


void SacRec::merge( SacRec sacrec2 ) {
   // make sure that both signals are loaded
   if( !sig || !sacrec2.sig )
		throw ErrorSR::EmptySig(FuncName, "for either sr1 or sr2");
   
   SAC_HD& shd2 = sacrec2.shd;
   // starting and ending time
   double t1b = this->AbsTime() + shd.b;
   double t2b = sacrec2.AbsTime() + sacrec2.shd.b;
   double t1e = t1b + (shd.npts-1)*shd.delta;
   double t2e = t2b + (shd2.npts-1)*shd2.delta;
   double T1 = std::min(t1b, t2b), T2 = std::max(t1e, t2e);

   double dt = (int)floor(shd.delta*1e8+0.5)/1e8, tshift;
   int N = (int)floor((T2-T1)/shd.delta+0.5)+1;

   /* pre-merging checks */
   if( N > 2.6e7 ) 
		throw ErrorSR::MemError( FuncName, "> 100M if merged");
   if( (shd.delta-shd2.delta) > 0.0001 )
		throw ErrorSR::BadParam( FuncName, "sps mismatch");

   /* allocate new space */
   std::unique_ptr<float[]> sig0(new float[N]);
	if( ! sig0 )
		throw ErrorSR::MemError( FuncName, "new failed!");
   std::fill(&(sig0[0]), &(sig0[N]), std::numeric_limits<float>::max()); // initialize the array to max float
   //for (j=0;j<N;j++) sig0[j] = 1.e30;
   //std::copy(&(sig1[0]), &(sig1[0])+5, &(sig0[0]));

   /* compute merge locations */
   int nb;
   bool reversed;
   if( t1b > t2b ) {
      reversed = true;
      nb = (int)floor((t1b-t2b)/dt+0.5);
		tshift = (nb*dt-(t1b-t2b));// + (shd.b-shd2.b);
   }
   else {
      reversed = false;
      nb = (int)floor((t2b-t1b)/dt+0.5);
		tshift = (nb*dt-(t2b-t1b));// + (shd2.b-shd.b);
   }
   if( fabs(tshift) > 1.e-3 )
		(*report) << "signal shifted by " << tshift << "secs" << std::endl;


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
	shd.e = shd.b + (N-1)*shd.delta;
}

int SacRec::arrange(const char* recname) {
   // count holes
   float maxfloat = std::numeric_limits<float>::max();
   int Nholes=0;
	float *sigsac = sig.get();
   for(int i=0;i<shd.npts;i++)
      if(sigsac[i]>=maxfloat) Nholes++;

   // produce record file for holes if a recname is given
   if( recname ) {
      int rec_b[1000], rec_e[1000];
      //char recname[200];
      rec_b[0] = 0;
      int j=0;
      for(int i=1;i<shd.npts;i++) {
         if(sigsac[i-1] >= maxfloat) { if(sigsac[i] < maxfloat) rec_b[j] = i; }
         else if(sigsac[i] >= maxfloat) rec_e[j++] = i;
      }
      if(sigsac[shd.npts-1]<maxfloat) rec_e[j++] = shd.npts;
      //sprintf(recname, "%s_rec1", ffsac.c_str());
      std::ofstream frec(recname);
      for(int i=0;i<j;i++) frec << rec_b[i] << " " << rec_e[i] << std::endl;
      //fprintf(frec, "%d %d\n", rec_b[i], rec_e[i]);
      frec.close();
   }

   // fill gaps with random numbers
   int ib, hlen = 100./shd.delta, step = std::max(0.5/shd.delta, 1.);
	bool isgap = false;
	for (int i=0; i<shd.npts; i++ ) {
		if( isgap ) {
			if ( sigsac[i]<maxfloat || i==shd.npts-1 ) {
				// compute rms in the surrounding time
				float sigtime_b = shd.b + ib * shd.delta;
				float sigtime_e = shd.b + i * shd.delta;
				float meanhlen = 500., sigmean1, sigmean2, sigstd;
				for( int ihlen=1; ihlen<(shd.e-shd.b)/meanhlen; ihlen++ ) {
					float rmswb = std::max( sigtime_b - meanhlen * ihlen, shd.b );
					float rmswe = std::min( sigtime_e + meanhlen * ihlen, shd.e );
					//if( RMSAvg( rmswb, rmswe, sigrms ) ) break;
					if( MeanStd( rmswb, rmswe, step, sigmean1, sigstd ) ) {
						float ftmp;
						sigmean1 = NaN; sigmean2 = NaN;
						MeanStd( rmswb, sigtime_b, step, sigmean1, ftmp );
						MeanStd( sigtime_e, rmswe, step, sigmean2, ftmp );
						if( sigmean1 == NaN ) sigmean1 = sigmean2;
						else if (sigmean2 == NaN ) sigmean2 = sigmean1;
						break;
					}
				}
				pimpl->FillGap(&(sig[ib]), &(sig[i]), sigmean1, sigmean2, sigstd, hlen, step );
				isgap = false;
			}
		}
		else {
			if ( sigsac[i] >= maxfloat ) {
				ib = i;
				isgap = true;
			}
		}
	}

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


/* ---------------------------------------- cut by event ---------------------------------------- */
void SacRec::ZoomToEvent( const std::string etime, float evlon, float evlat, float tb, float tlen, std::string ename ) {
	if( etime.length() != 14 )
		throw ErrorSR::BadParam( FuncName, "Invalid etime string = " + etime );
	SAC_HD eshd;
   eshd.nzyear = std::stoi( etime.substr(0, 4) );
   float month = std::stoi( etime.substr(4, 2) );
   float day = std::stoi( etime.substr(6, 2) );
   eshd.nzjday = pimpl->Jday( eshd.nzyear, month, day );
   eshd.nzhour = std::stoi( etime.substr(8, 2) );
   eshd.nzmin = std::stoi( etime.substr(10, 2) );
   eshd.nzsec = std::stoi( etime.substr(12, 2) );
   eshd.nzmsec = 0.;
	if( ename.empty() ) ename = etime;
	ZoomToEvent( eshd, evlon, evlat, tb, tlen, ename );
}

void SacRec::ZoomToEvent( const SAC_HD& eshd, float evlon, float evlat, float tb, float tlen, const std::string ename ) {
   if( tlen <= 0. )
		throw ErrorSR::BadParam( FuncName, "negative tlen = " + std::to_string(tlen) );
   
   // current origin time
	double Tabs_old = AbsTime();

   // assign event time to sacrec and update shd.b
   shd.nzyear = eshd.nzyear; shd.nzjday = eshd.nzjday;
   shd.nzhour = eshd.nzhour; shd.nzmin = eshd.nzmin;
	shd.nzsec = eshd.nzsec; shd.nzmsec = eshd.nzmsec;
	double Tabs_new = AbsTime();
   shd.b -= ( Tabs_new - Tabs_old );
	shd.e = shd.b + (shd.npts-1) * shd.delta;

   // cut to tb - te
	float te = tb + tlen;
   cut( tb, te );

   // assign event location
   sprintf( shd.kevnm, "%s", ename.c_str() );
	if( evlon>=-180. && evlon<=360. ) shd.evlo = evlon;
   if( evlat>=-90. && evlat<=90. ) shd.evla = evlat;
}

/* remove mean and trend */
void SacRec::RTrend() {
   // fit a*x+b
   int i, npts = shd.npts;
   float X = 0., Y = 0., X2 = 0., Y2 = 0., XY = 0.;
	float *sigsac = sig.get();
   for(i=0;i<npts;i++) {
      X += i;
      Y += sigsac[i];
      X2 += i*i;
      Y2 += sigsac[i]*sigsac[i];
      XY += i*sigsac[i];
   }
   float a = (npts*XY-X*Y)/(npts*X2-X*X);
   float b = (-X*XY+X2*Y)/(npts*X2-X*X);
   // correct sig and DEPMEN
   float mean = 0., max = sigsac[0], min = sigsac[0];
   float shift = b;
   for(i=0;i<npts;i++,shift+=a) {
      sigsac[i] -= shift;
      mean += sigsac[i];
      if ( min > sigsac[i] ) min = sigsac[i];
      else if ( max < sigsac[i] ) max = sigsac[i];
   }
   shd.depmin = min;
   shd.depmax = max;
   shd.depmen = mean / npts;
}


/* remove response and apply filter */
void SacRec::RmRESP( const std::string& fresp, float perl, float perh, const std::string& evrexein, const int type ) {
	
   // check evrexe
	std::string evrexe(evrexein);
   if( evrexe.empty() )
      pimpl->FindInPath("evalresp", evrexe);
   if( access( evrexe.c_str(), F_OK ) == -1 ) 
		throw ErrorSR::BadParam( FuncName, "cannot access evralresp at " + evrexe );

   // run evalresp
   int nf = 100;
   char buff[300], sta[8], ch[8], net[8];
   float f2 = 1./perh, f1 = f2*0.9, f3 = 1./perl, f4 = f3*1.1;
   //float f2 = 1./perh, f1 = f2*0.8, f3 = 1./perl, f4 = f3*1.2;
   sscanf(shd.kstnm, "%s", sta);
   sscanf(shd.kcmpnm, "%s", ch);
   sscanf(shd.knetwk, "%s", net);
   sprintf(buff, "%s %s %s %4d %3d %f %f %d -f %s -v >& /dev/null", evrexe.c_str(), sta, ch, shd.nzyear, shd.nzjday, f1, f4, nf, fresp.c_str());
	// define all variables to be used in the critical section
   char nameam[50], nameph[50];
   sprintf(nameam, "AMP.%s.%s.*.%s", net, sta, ch);
   sprintf(nameph,"PHASE.%s.%s.*.%s", net, sta, ch);
   double freq[nf], amp[nf], pha[nf];
	bool openA=false, openP=false;
	int sizeA=0, sizeP=0;
	#pragma omp critical(external)
	{
   system(buff);
   // find am file
   FILE *fam = nullptr, *fph = nullptr;
   std::vector<std::string> list;
   System::List(".", nameam, 0, list);
	sizeA = list.size();
	if( sizeA == 1 ) {
		sscanf(list.at(0).c_str(), "%s", nameam);
		fam = fopen(nameam, "r");
		if( fam ) openA = true;
	}
   // find ph file
   System::List(".", nameph, 0, list);
	sizeP = list.size();
   if( sizeP == 1 ) {
		sscanf(list.at(0).c_str(), "%s", nameph);
		fph = fopen(nameph, "r");
		if( fph ) openP = true;
	}
	// proceed further only if all succed!
	if( sizeA==1 && sizeP==1 && openA && openP ) {
		// read in am and ph data
		double pi=4*atan(1.0), pio180=pi/180., dtmp;
		int i = 0;
		while(i<nf) {
			if(fgets(buff, 300, fam)==nullptr) break;
			sscanf(buff, "%lf %lf", &freq[i], &amp[i]);
			if(fgets(buff, 300, fph)==nullptr) break;
			sscanf(buff, "%lf %lf", &dtmp, &pha[i]);
			if(dtmp!=freq[i]) {
				(*report)<<"incompatible AMP - PHASE pair!"<<std::endl;
				continue;
			}
			amp[i] *= 0.000000001;
			pha[i] *= pio180;
			i++;
		}
	}
	if( openA ) fclose(fam); 
	if( openP ) fclose(fph);
	if( sizeA==1 ) System::fRemove(nameam); 
	if( sizeP==1 ) System::fRemove(nameph);
	}
	// throw exceptions if failed
	if( sizeA != 1 )
		throw ErrorSR::ExternalError( FuncName, std::to_string(sizeA) + " AMP file(s) found");
	if( ! openA )
		throw ErrorSR::BadFile( FuncName, "reading from " + std::string(nameam) );
	if( sizeP != 1 )
		throw ErrorSR::ExternalError( FuncName, std::to_string(sizeP) + " PHASE file(s) found");
	if( ! openP )
		throw ErrorSR::BadFile( FuncName, "reading from " + std::string(nameph) );
   // remove trend ( and mean )
   RTrend();
   // run rmresponse
	if( shd.npts > maxnpts4parallel ) {
		#pragma omp critical(largesig)
		{
		pimpl->FDivide (f1, f2, f3, f4, freq, amp, pha, nf, shd, sig.get());
		}
	} else {
		pimpl->FDivide (f1, f2, f3, f4, freq, amp, pha, nf, shd, sig.get());
	}

	// convert to the requested data type:
	switch(type) {
		case 0:
			Integrate(); break;
		case 1:
			break;
		case 2:
			Differentiate(); break;
		default:
			throw ErrorSR::BadParam( FuncName, "unknown datatype (type="+std::to_string(type)+") for sac output." );
	}
}

/* down sampling with anti-aliasing filter */
void SacRec::Resample( int sps, bool fitParabola ) {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);
	if( sps <= 0 ) sps = floor(1.0/shd.delta+0.5);
		
	// grid step size
   float dt = 1./sps;
	if( dt < shd.delta )
		throw ErrorSR::BadParam( FuncName, "Upsampling not implemented" );
   int iinc = (int)floor(dt/shd.delta+0.5);
	
   // anti-aliasing filter
   if(iinc!=1) {
      double f3 = sps/2.2, f4 = sps/2.01;
		LowpassCOSFilt(f3, f4);
   }

	// calculate starting grid
   float t0 = shd.nzmsec*0.001+shd.b, nb = ceil(t0*sps), t0new = nb * dt;
   int i = (int)floor((t0new-t0)/shd.delta);
	if( iinc==1 && t0==0. ) return;	// nothing needs to be done

   // allocate space for the new sig pointer
	int nptst = (int)floor( ( t0 + (shd.npts-1)*shd.delta - t0new ) * sps ) + 1;
   //int nptst = nint((shd.npts-1)*shd.delta*sps)+10;
   std::unique_ptr<float[]> sig2(new float[nptst]);
	if( ! sig2 )
		throw ErrorSR::MemError( FuncName, "new failed!");
	float *sigsac = sig.get(), *sigsac2 = sig2.get();

	// re-sample
   //if( (*sig2 = (float *) malloc (nptst * sizeof(float)))==nullptr ) perror("malloc sig2");
   int j;
	if( fitParabola ) {	// slower yet more accurate
		shd.b += shd.nzmsec*0.001;
      //long double ti = i*shd.delta + t0;
      long double tj = t0new;
		for(j=0; j<nptst; j++) {
			sigsac2[j] = Sig(tj);
			tj += dt;
		}
   } else if(fabs(iinc*shd.delta-dt)<1.e-7) { //sps is a factor of 1/delta
      long double fra2 = (t0new-t0-i*shd.delta)/shd.delta;
      long double fra1 = 1.-fra2;
      if(fra2==0)
         for(j=0;i<shd.npts;j++) {
            sigsac2[j] = sigsac[i];
            i += iinc;
         }
      else
         for(j=0;i<shd.npts-1;j++) {
            sigsac2[j] = sigsac[i]*fra1 + sigsac[i+1]*fra2;
            i += iinc;
         }
   } else { //sps isn't a factor, slower way
      (*report) << "possible rounding error! sps isn't a factor of " << (int)floor(1/shd.delta+0.5) << std::endl;
      long double ti, tj;
      iinc = (int)floor(dt/shd.delta);
      ti = i*shd.delta + t0;
      tj = t0new;
      for(j=0;i<shd.npts-1;j++) {
         long double fra2 = tj-ti;
         sigsac2[j] = sigsac[i] + (sigsac[i+1]-sigsac[i])*fra2;
         tj += dt;
         i += iinc;
         ti += iinc*shd.delta;
         if( ti+shd.delta <= tj ) { ti += shd.delta; i++; }//if(j%1000==0)cerr<<i<<" "<<ti<<j<<" "<<tj<<" "<<endl;}
      }
   }
   sig = std::move(sig2);
   //shd.b = 0.;
   //shd.nzmsec = (int)(nb*dt*1000+0.5);
   //if(shd.nzmsec>=1000) UpdateTime();
	// modify shd.b instead of the KZtime
	shd.nzmsec = 0.;
	shd.b = t0new;
   shd.delta = dt;
   shd.npts = j;
}



/* ---------------------------------------- temporal normalizations ---------------------------------------- */
void SacRec::OneBit() {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);

	float *sigsac = sig.get();
   for(int i=0;i<shd.npts;i++) {
      if(sigsac[i]>0.) sigsac[i] = 1.;
      else if(sigsac[i]<0.) sigsac[i] = -1.;
   }
}


void SacRec::RunAvg( float timehlen, float Eperl, float Eperh ) {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);

	/* filter into the earthquake band */
	SacRec sac_eqk;
	if( Eperl == -1. ) {
		sac_eqk = *this;
	} else {
		float f2 = 1./Eperh, f1 = f2*0.6, f3 = 1./Eperl, f4 = f3*1.4;
		BandpassCOSFilt( f1, f2, f3, f4, sac_eqk );
	}

	/* smooth to get earthquake strength */
   int i, j, wb, we, n = shd.npts;
   float wsum, dt = shd.delta;
   int half_l = (int)floor(timehlen/dt+0.5);
	float *sigsac = sig.get(), *sigw = sac_eqk.sig.get();
	for( i=0; i<shd.npts; i++ ) sigw[i] = fabs(sigw[i]);

   if(half_l*2>n-1) half_l = (n-1)/2;
   for(i=0,wsum=0.;i<=half_l;i++) wsum += sigw[i];
   wb = 0; we = i;
   for(i=1;i<=half_l;i++,we++) {
      if(wsum>1.e-15) sigsac[i-1] *= ((double)we/wsum);
      wsum += sigw[we];
   }
   for(j=we;i<n-half_l;i++,wb++,we++) {
		//if( i>80000/dt && i<82000/dt ) std::cerr<<(i-1)*dt<<" "<<sig[i-1]<<" "<<wsum<<" / "<<j<<std::endl;
      if(wsum>1.e-15) sigsac[i-1] *= ((double)j/wsum); //fout<<shd.b+(i-1)*dt<<" "<<wsum<<" "<<j<<std::endl;}
      wsum += ( sigw[we] - sigw[wb] );
   }
   for(;i<n;i++,wb++) {
      if(wsum>1.e-15) sigsac[i-1] *= ((double)(we-wb)/wsum);
      wsum -= sigw[wb];
   }
   if(wsum>1.e-15) sigsac[n-1] *= ((double)(we-wb)/wsum);


}

/* ---------- compute the correlation coefficient with an input SacRec ---------- */
float SacRec::Correlation( const SacRec& sac2, const float tb, const float te ) const {
	// check input
	if( !sig || !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( shd.delta != sac2.shd.delta )
		throw ErrorSR::HeaderMismatch(FuncName, "dt");

	const auto& sac1 = *this;

	// compute mean and std-devs
	bool succeed = true;
	float mean1, mean2, std1, std2;
	succeed &= sac1.MeanStd( tb, te, mean1, std1 );
	succeed &= sac2.MeanStd( tb, te, mean2, std2 );
	if( ! succeed )
		throw ErrorSR::BadParam(FuncName, "time window out of range: " + std::to_string(tb)+" - "+std::to_string(te));

	// get/check indexes
	float *sigsac1 = sac1.sig.get();
	float *sigsac2 = sac2.sig.get();
	float cc = 0.;
	size_t ib1 = Index(tb), ie1 = Index(te) + 1;
	size_t ib2 = sac2.Index(tb), ie2 = sac2.Index(te) + 1;
	if( X(ib1) != sac2.X(ib2) )
		throw ErrorSR::HeaderMismatch(FuncName, "sampling grid");

	// compute correlation coef
	size_t i1 = ib1, i2 = ib2;
	while( i1 < ie1 ) {
		cc += (sigsac1[i1]-mean1) * (sigsac2[i2]-mean2);
		i1++; i2++;
	}
	//std::cerr<<cc<<"   "<<tb<<" "<<ib<<" "<<te<<" "<<ie<<"   "<<mean1<<" "<<std1<<" "<<mean2<<" "<<std2<<std::endl;
	cc /= ( (ie1-ib1-1) * std1 * std2 );

	return cc;
}

/* Cross-Correlate with another sac record
	ctype=0: Cross-Correlate (default) 
	ctype=1: deconvolve (sac.am/sac2.am)
	ctype=2: deconvolve (sac2.am/sac.am) */
/*
void SacRec::CrossCorrelate( SacRec& sac2, SacRec& sacout, int ctype ) {
	// check input
	if( !sig || !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::SizeMismatch(FuncName, std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );

	// sig1 time -> freq
	SacRec sac1_am, sac1_ph;
	ToAmPh(sac1_am, sac1_ph);
	
	// sig2 time -> freq
	SacRec sac2_am, sac2_ph;
	sac2.ToAmPh(sac2_am, sac2_ph);

	CrossCorrelate( sac1_am, sac1_ph, sac2_am, sac2_ph, sacout, ctype );

}
*/
void SacRec::CrossCorrelate( SacRec& sac2, SacRec& sacout, int ctype ) {
	// check input
	if( !sig || !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( shd.delta != sac2.shd.delta )
		throw ErrorSR::HeaderMismatch(FuncName, "dt");
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );

	// sig1 time -> freq
	SacRec sac1_am, sac1_ph;
	ToAmPh(sac1_am, sac1_ph);
	
	// sig2 time -> freq
	SacRec sac2_am, sac2_ph;
	sac2.ToAmPh(sac2_am, sac2_ph);

	// initialize amp and phase out
	int npts = sac1_am.shd.npts;
	SacRec out_am;
	if( ctype == 2 ) out_am = sac2_am;
	else out_am = sac1_am;
	/*
	out_am.sig.reset( new float[ns]() );
	if( ! out_am.sig )
		throw ErrorSR::MemError( FuncName, "new failed for out_am!");
	out_am.shd = sac1_am.shd;
	out_am.shd.npts = ns;
	*/

	// correlate in freq domain
	float *amout = out_am.sig.get();
	float *am1 = sac1_am.sig.get(), *am2 = sac2_am.sig.get();
	float ampmin = 0.;
	switch( ctype ) {
		case 0:
			for(int i=0; i<npts; i++) amout[i] *= am2[i];
			break;
		case 1:
			for(int i=0; i<npts; i++) ampmin += am2[i];
			ampmin /= npts*20.;
			for(int i=0; i<npts; i++) amout[i] /= std::max(am2[i], ampmin);
			break;
		case 2:
			for(int i=0; i<npts; i++) ampmin += am1[i];
			ampmin /= npts*20.;
			for(int i=0; i<npts; i++) amout[i] /= std::max(am1[i], ampmin);
			break;
		default:
			throw ErrorSR::BadParam( FuncName, "Unknown Cor type!" );
	}

	// clear sacs & free memories 1
	sac1_am.clear(); sac2_am.clear();

	// compute phase out
	SacRec out_ph( sac2_ph );
	/*
	out_ph.sig.reset( new float[ns]() );
	if( ! out_ph.sig )
		throw ErrorSR::MemError( FuncName, "new failed for out_ph!");
	out_ph.shd = sac1_ph.shd;
	out_ph.shd.npts = ns;
	*/
	float *phout = out_ph.sig.get();
	float *ph1 = sac1_ph.sig.get(), *ph2 = sac2_ph.sig.get();
   for(int i=0; i<npts; i++) phout[i] -= ph1[i];

	// clear sacs & free memories 2
	sac1_ph.clear(); sac2_ph.clear();

	// convert back to time domain (with length ns)
	SacRec sac_CC;
	sac_CC.FromAmPh(out_am, out_ph);
	// clear sacs & free memories 3
	out_am.clear(); out_ph.clear();

	// form final result
	int lag = shd.npts-1, ns = (npts-1) * 2;
	sacout.shd = shd;
	sacout.shd.npts = lag*2 + 1;
	float Tshift = sac2.shd.b-shd.b, Tlag = lag * sacout.shd.delta;
	sacout.shd.b = -Tlag + Tshift;
	sacout.shd.e = Tlag + Tshift;
	sacout.sig.reset( new float[lag*2+1]() );
	float *cor = sacout.sig.get(), *CCout = sac_CC.sig.get();
   for( int i = 1; i< (lag+1); i++) {
      cor[lag+i] =  CCout[i];
      cor[lag-i] =  CCout[ns-i];
   }
   cor[lag] = CCout[0];

	// normalize
	//float *outsig = sacout.sig.get();
   //for(int i=0; i<sacout.npts; i++) outsig[i] *= ;
	
}


