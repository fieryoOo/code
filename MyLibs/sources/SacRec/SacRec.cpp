#include "SacRec.h"
#include "MyLogger.h"
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
   #define PI 3.14159265358979323846
   // forward FFT 
   void FFTW_F(fftw_plan plan, fftw_complex *out, int ns, float *seis, int n) {
      fftw_execute(plan);
      //pthread_mutex_lock(&fftlock);
		#pragma omp critical(fftw)
		{
      fftw_destroy_plan(plan);
		}
      //pthread_mutex_unlock(&fftlock);
      int k;
      for(k=0; k<n; k++) seis[k] = out[k][0];
      //for(k=0; k<n; k+=1000) if(seis[k] != 0.) printf("%d %f\n", k, seis[k]);
   }
   // backward FFT
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
      TaperB( f1, f2, f3, f4, dom, nk, sf );

      //forward FFT: sf ==> s
      //fftw_plan plan2;
      FFTW_F(plan1, s, ns, seis_in, n);
      fftw_free(s); fftw_free(sf);

      //forming final result
      int k;
      float ftmp = 2./ns;
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


extern MyLogger logger;
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
   int npts=recin.shd.npts; sig.reset(new float[npts]);
	if( ! sig )
		throw ErrorSR::MemError( FuncName, "new failed!");
   std::copy(recin.sig.get(), recin.sig.get()+npts, sig.get());
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
   fsac.read( reinterpret_cast<char *>(&shd), sizeof(SAC_HD) );
   fsac.close();
   //pthread_mutex_unlock(&fiolock);
}

/* load sac header+signal from file 'fname', memory is allocated on heap */
void SacRec::Load () {
   std::ifstream fsac(fname.c_str());
   if( ! fsac )
		throw ErrorSR::BadFile( FuncName, "reading from " + fname );
   //if( SHDMap.empty() ) pimpl->CreateSHDMap();
   //pthread_mutex_lock(&fiolock);
   fsac.read( reinterpret_cast<char *>(&shd), sizeof(SAC_HD) );
   //sig = std::make_shared<float>( new float[shd.npts*sizeof(float)] );
   float* sigtmp = new float[shd.npts];
	if( sigtmp == nullptr )
		throw ErrorSR::MemError( FuncName, "new failed!");
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
}

/* write to file '*outfname' */
void SacRec::WriteHD (const std::string& outfname) {
   /* open file */
   //std::fstream fsac(outfname, std::ios::in | std::ios::out);
   std::fstream fsac(outfname);
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

   /* search min and max amplitude */
	float *sigsac = sig.get();
   shd.depmin = sigsac[0];
   shd.depmax = sigsac[0];
   for ( int i = 0; i < shd.npts ; i++ ) {
       if ( shd.depmin > sigsac[i] ) shd.depmin = sigsac[i];
       else if ( shd.depmax < sigsac[i] ) shd.depmax = sigsac[i];
   }

   /* check and re-format header time if necessary */
   UpdateTime();

   //pthread_mutex_lock(&fiolock);
   fsac.write( reinterpret_cast<char *>(&shd), sizeof(SAC_HD) );
   fsac.write( reinterpret_cast<char *>(sigsac), sizeof(float)*shd.npts );
   //pthread_mutex_unlock(&fiolock);

   fsac.close();
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
//#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
//#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

//pthread_mutex_t fiolock;

void SacRec::Mul( const float mul ) {
	float *sigsac = sig.get();
   for(int i=0; i<shd.npts; i++) sigsac[i] *= mul;
}

void SacRec::Addf( const SacRec& sac2 ) {
   if( !sig || !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::SizeMismatch(FuncName, std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );
	//const auto& sig2 = sac2.sig;
	float *sigsac = sig.get(), *sigsac2 = sac2.sig.get();
   for(int i=0; i<shd.npts; i++) 
		sigsac[i] += sigsac2[i];
}

void SacRec::Divf( const SacRec& sac2 ) {
   if( !sig || !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::SizeMismatch(FuncName, std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );
	//const auto& sig2 = sac2.sig;
	float *sigsac = sig.get(), *sigsac2 = sac2.sig.get();
   for(int i=0; i<shd.npts; i++) 
		if( sigsac2[i]!=0. ) sigsac[i] /= sigsac2[i];
}

void SacRec::PullUpTo( const SacRec& sac2 ) {
   if( !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( !sig )
		*this = sac2;
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::SizeMismatch(FuncName, std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );
	//const auto& sig2 = sac2.sig;
	float *sigsac = sig.get(), *sigsac2 = sac2.sig.get();
   for(int i=0; i<shd.npts; i++) 
		if( sigsac2[i] > sigsac[i] ) sigsac[i] = sigsac2[i];
		//sigsac[i] = std::max(sigsac[i], sigsac2[i]);
}


#include <cctype>
void SacRec::ChHdr(const std::string& fieldin, const std::string& value){
   std::stringstream sin(value);
	std::string field(fieldin);
   std::transform(field.begin(), field.end(), field.begin(), ::tolower);
   bool succeed = false;

   if( field == "dist" ) succeed = sin >> shd.dist;
   else if( field == "az" ) succeed = sin >> shd.az;
   else if( field == "baz" ) succeed = sin >> shd.baz;
   else if( field == "gcarc" ) succeed = sin >> shd.gcarc;
   else if( field == "b" ) succeed = sin >> shd.b;
   else if( field == "e" ) succeed = sin >> shd.e;

   else if( field == "knetwk" ) succeed = sin >> shd.knetwk;
   else if( field == "kstnm" ) succeed = sin >> shd.kstnm;
   else if( field == "stlo" ) succeed = sin >> shd.stlo;
   else if( field == "stla" ) succeed = sin >> shd.stla;
   else if( field == "stel" ) succeed = sin >> shd.stel;
   else if( field == "stdp" ) succeed = sin >> shd.stdp;

   else if( field == "kevnm" ) succeed = sin >> shd.kevnm;
   else if( field == "evlo" ) succeed = sin >> shd.evlo;
   else if( field == "evla" ) succeed = sin >> shd.evla;
   else if( field == "evel" ) succeed = sin >> shd.evel;
   else if( field == "evdp" ) succeed = sin >> shd.evdp;

   else if( field == "nzyear" ) succeed = sin >> shd.nzyear;
   else if( field == "nzjday" ) succeed = sin >> shd.nzjday;
   else if( field == "nzhour" ) succeed = sin >> shd.nzhour;
   else if( field == "nzmin" ) succeed = sin >> shd.nzmin;
   else if( field == "nzsec" ) succeed = sin >> shd.nzsec;
   else if( field == "nzmsec" ) succeed = sin >> shd.nzmsec;

   else if( field == "kcmpnm" ) succeed = sin >> shd.kcmpnm;
   else if( field == "cmpaz" ) succeed = sin >> shd.cmpaz;
   else if( field == "cmpinc" ) succeed = sin >> shd.cmpinc;

   else if( field == "o" ) succeed = sin >> shd.o;
   else if( field == "ko" ) succeed = sin >> shd.ko;
   else if( field == "a" ) succeed = sin >> shd.a;
   else if( field == "ka" ) succeed = sin >> shd.ka;
   else if( field == "f" ) succeed = sin >> shd.f;
   else if( field == "kf" ) succeed = sin >> shd.kf;

   else if( field == "user0" ) succeed = sin >> shd.user0;
   else if( field == "user1" ) succeed = sin >> shd.user1;
   else if( field == "user2" ) succeed = sin >> shd.user2;
   else if( field == "user3" ) succeed = sin >> shd.user3;
   else if( field == "user4" ) succeed = sin >> shd.user4;
   else if( field == "user5" ) succeed = sin >> shd.user5;
   else if( field == "user6" ) succeed = sin >> shd.user6;
   else if( field == "user7" ) succeed = sin >> shd.user7;
   else if( field == "user8" ) succeed = sin >> shd.user8;
   else if( field == "user9" ) succeed = sin >> shd.user9;

	if( ! succeed )
		throw ErrorSR::BadParam(FuncName, "Invalid header field (" + field + ") or value (" + value + ")");
}


double SacRec::AbsTime () {
   //if( ! sig ) return -1.;
   //if( shd == sac_null ) return -1.; // operator== not defined yet
   if( shd.npts <= 0 ) return -1.;
   if( shd.nzjday == -12345. || shd.nzyear == -12345. || shd.nzhour == -12345. ||
       shd.nzmin == -12345. || shd.nzsec == -12345. || shd.nzmsec == -12345. ) return -1;
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
void SacRec::MinMax (float tbegin, float tend, float& tmin, float& min, float& tmax, float& max) {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);
	float *sigsac = sig.get();
   min = max = sigsac[0];
   int imin = 0, imax = 0;
   for ( int i = nint((tbegin-shd.b)/shd.delta); i < nint((tend-shd.b)/shd.delta+1.) ; i++ ) {
       if ( min > sigsac[i] ) { min = sigsac[i]; imin = i; }
       else if ( max < sigsac[i] ) { max = sigsac[i]; imax = i; }
   }
   tmin = shd.b + imin*shd.delta;
   tmax = shd.b + imax*shd.delta;
}


/* compute the root-mean-square average in a given window */
void SacRec::RMSAvg ( float tbegin, float tend, int step, float& rms ) {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);
	/* decide avg window */
   float maxfloat = std::numeric_limits<float>::max();
   int neff = 0, ibeg = nint((tbegin-shd.b)/shd.delta), iend = nint((tend-shd.b)/shd.delta+1.);
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

   rms = 0.;
	float *sigsac = sig.get();
   for ( int i = ibeg; i < iend ; i+=step ) {
      if( sigsac[i] >= maxfloat ) continue;
      rms += sigsac[i] * sigsac[i];
      neff++;
   }
   rms = sqrt(rms/(neff-1));
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
bool SacRec::MeanStd ( float tbegin, float tend, int step, float& mean, float& std ) {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);

	// store valid data points
   float maxfloat = std::numeric_limits<float>::max();
   int neff = 0, ibeg = nint((tbegin-shd.b)/shd.delta), iend = nint((tend-shd.b)/shd.delta+1.);
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


/* ---------------------------------------- time - frequency ---------------------------------------- */
/* convert to amplitude */
void SacRec::ToAmPh( SacRec& sac_am, SacRec& sac_ph ) {
	if( shd.npts > maxnpts4parallel ) {
		#pragma omp critical(largesig)
		{
		ToAmPh_p(sac_am, sac_ph);
		}
	} else {
		ToAmPh_p(sac_am, sac_ph);
	}
}
void SacRec::ToAmPh_p( SacRec& sac_am, SacRec& sac_ph ) {
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
   int nk = ns/2 + 1;
   float *amp = new float[nk], *pha = new float[nk];
	if( amp==nullptr || pha==nullptr )
		throw ErrorSR::MemError( FuncName, "new failed!");
   for(int i=0; i<nk; i++) {
      fftw_complex& cur = sf[i];
      amp[i] = sqrt(cur[0]*cur[0] + cur[1]*cur[1]);
		pha[i] = atan2(cur[1], cur[0]);
   }
   //fftw_free(s); 
	fftw_free(sf);
   sac_am.sig.reset(amp);
	sac_ph.sig.reset(pha);

   sac_am.shd.npts = nk;
   sac_am.shd.delta = 1./(shd.delta*ns);
   sac_am.shd.b = 0.;
	sac_ph.shd = sac_am.shd;
}

/* 3 different types of filters */
void SacRec::Filter ( double f1, double f2, double f3, double f4, SacRec& srout ) {
	if( shd.npts > maxnpts4parallel ) {
		#pragma omp critical(largesig)
		{
		Filter_p(f1, f2, f3, f4, srout);
		}
	} else {
		Filter_p(f1, f2, f3, f4, srout);
	}
}
void SacRec::Filter_p ( double f1, double f2, double f3, double f4, SacRec& srout ) {
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
   if(f4 > 0.5/dt) {
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
   int nk = ns/2+1;
   double dom = 1./(dt*ns);
   if( (f1==-1. || f2==-1.) && (f3>0. && f4>0.) ) pimpl->TaperL( f3, f4, dom, nk, sf );
   else if( f1>=0. && f2>0. && f3>0. && f4>0. ) pimpl->TaperB( f1, f2, f3, f4, dom, nk, sf );
   else if( f1==-1. && f4==-1. ) pimpl->TaperGaussian( f2, f3, dom, nk, sf );
   else throw ErrorSR::BadParam( FuncName, "Unknown filter type");

   //forward FFT: sf ==> s
   pimpl->FFTW_F(planF, s, ns, &(srout.sig[0]), n);
   fftw_free(s); fftw_free(sf);

   //forming final result
   float ftmp = 2./ns;
	float *sigout = srout.sig.get();
   for(int k=0; k<n; k++) sigout[k] *= ftmp;
   //   if( seis_in[k]==0 ) seis_out[k] = 0.;
   //   else seis_out[k] *= ftmp;
}


/* cosine tapper */
void SacRec::cosTaperL( const float fl, const float fh ) {
   if( !sig || shd.npts<=0 )	// check signal
		throw ErrorSR::EmptySig(FuncName);

	const float& delta = shd.delta;
   int i;
	float *sigsac = sig.get();
	for( i=0; i<(int)ceil(fl/delta); i++ ) sigsac[i] = 0.;
	float finit = i*delta, fwidth = fh-fl, dfh = fh-finit;
	float ftmp = PI / fwidth;
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
   int i = (int)ceil(fl/delta);
	float finit = i*delta, fwidth = fh-fl, dfl = finit-fl;
	float ftmp = PI / fwidth;
	float *sigsac = sig.get();
   for(; dfl<fwidth; i++, dfl+=delta) {
		float amplif = ( 1. + cos(ftmp*dfl) ) * 0.5;
		sigsac[i] *= amplif;
   }
   for(;i<shd.npts;i++) sigsac[i] = 0.;

   return;
}


/* ---------------------------------------- cut and merge ---------------------------------------- */
void SacRec::cut( float tb, float te, SacRec& sac_result ) {
   int nb = (int)floor( (tb-shd.b) / shd.delta + 0.5 );
   int ne = (int)floor( (te-shd.b) / shd.delta + 0.5 );
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
   float nptscpy = std::min( nptsnew - inew, shd.npts - iold) - 1;
   // copy data
   memcpy( &(signew[inew]), &(sig[iold]), nptscpy * sizeof(float) );
   // reset sacT.sig
   sac_result.sig.reset(signew);
   // update shd
   if( &(sac_result) != this ) sac_result.shd = shd;
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
   double t1b = this->AbsTime();
   double t2b = sacrec2.AbsTime();
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
      tshift = (shd.b-shd2.b) + (nb*dt-(t1b-t2b));
   }
   else {
      reversed = false;
      nb = (int)floor((t2b-t1b)/dt+0.5);
      tshift = (shd2.b-shd.b) + (nb*dt-(t2b-t1b));
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
						sigmean1 = -12345.; sigmean2 = -12345.;
						MeanStd( rmswb, sigtime_b, step, sigmean1, ftmp );
						MeanStd( sigtime_e, rmswe, step, sigmean2, ftmp );
						if( sigmean1 == -12345. ) sigmean1 = sigmean2;
						else if (sigmean2 == -12345. ) sigmean2 = sigmean1;
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
void SacRec::RmRESP( const std::string& fresp, float perl, float perh, const std::string& evrexein ) {
	
   // check evrexe
	std::string evrexe(evrexein);
   if( evrexe.empty() )
      pimpl->FindInPath("evalresp", evrexe);
   if( access( evrexe.c_str(), F_OK ) == -1 ) 
		throw ErrorSR::BadParam( FuncName, "cannot access evralresp at " + evrexe );

   // run evalresp
   int nf = 100;
   char buff[300], sta[8], ch[8], net[8];
   //float f2 = 1./perh, f1 = f2*0.9, f3 = 1./perl, f4 = f3*1.1;
   float f2 = 1./perh, f1 = f2*0.8, f3 = 1./perl, f4 = f3*1.2;
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
}

/* down sampling with anti-aliasing filter */
void SacRec::Resample( float sps ) {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);
		
   int ithread = 0;
   /* anti-aliasing filter */
   float dt = 1./sps;
   int iinc = (int)floor(dt/shd.delta+0.5);
   if(iinc!=1) {
      double f3 = sps/2.2, f4 = sps/2.01;
		Filter(-1., 0., f3, f4);
   }

   /* allocate space for the new sig pointer */
   int i, j;
   int nptst = (int)floor((shd.npts-1)*shd.delta*sps+0.5)+10;
   float nb;
   std::unique_ptr<float[]> sig2(new float[nptst]);
	if( ! sig2 )
		throw ErrorSR::MemError( FuncName, "new failed!");
   //if( (*sig2 = (float *) malloc (nptst * sizeof(float)))==nullptr ) perror("malloc sig2");
   long double fra1, fra2;
	float *sigsac = sig.get(), *sigsac2 = sig2.get();
   nb = ceil((shd.nzmsec*0.001+shd.b)*sps);
   i = (int)floor((nb*dt-shd.nzmsec*0.001-shd.b)/shd.delta);
   if(fabs(iinc*shd.delta-dt)<1.e-7) { //sps is a factor of 1/delta
      fra2 = (nb*dt-i*shd.delta-shd.nzmsec*0.001-shd.b)/shd.delta;
      fra1 = 1.-fra2;
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
   }
   else { //sps isn't a factor, slower way
      (*report) << "possible rounding error! sps isn't a factor of " << (int)floor(1/shd.delta+0.5) << std::endl;
      long double ti, tj;
      iinc = (int)floor(dt/shd.delta);
      ti = i*shd.delta+shd.nzmsec*0.001+shd.b;
      tj = nb*dt;
      for(j=0;i<shd.npts-1;j++) {
         fra2 = tj-ti;
         sigsac2[j] = sigsac[i] + (sigsac[i+1]-sigsac[i])*fra2;
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
		Filter( f1, f2, f3, f4, sac_eqk );
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

