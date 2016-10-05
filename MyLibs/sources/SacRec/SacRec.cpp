#include "SacRec.h"
#include "DisAzi.h"
#include "Array2D.h"
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
		/* comparing plan to null. Is it a defined behavior?
      if( plan==nullptr || (planF && *planF==nullptr) ) {
         fprintf(stderr,"Error(FFTW_B): fftw_plan creation failed!!\n");
         //pthread_mutex_unlock(&fftlock); exit(0);
      } */
      //pthread_mutex_unlock(&fftlock);
      //initialize input array and excute
      memset(*in, 0, ns*sizeof(fftw_complex));
      int k;
      for(k=0; k<n; k++) (*in)[k][0] = seis[k];
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
      (*out)[0][0] *= 0.5; (*out)[0][1] *= 0.5;
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
			double s = 1. / std::sqrt(1. + pow(f/fc, index)), ss = s;
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
			double s = 1. / std::sqrt(1. + pow(fc/f, index)), ss = s;
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
	inline double btwB( double fcL, double fcR, int norder, double f ) {
		int index = 2 * norder;
		double s = 1. / std::sqrt(1. + pow((f*f-fcR*fcL)/((fcR-fcL)*f), index));
		return s;
	}
   void btwTaperB( double fcL, double fcR, int norder, double dom, int nk, fftw_complex *sf, int norderT ) {
		//btwTaperL( fcL, nL, dom, nk, sf );	// wrong
		//btwTaperR( fcR, nR, dom, nk, sf );	// and wrong!
		// apply taper
		int i; double f;
		int index = 2 * norder;
		for( i=0, f=0.; i<nk; i++, f+=dom) {
			double s = btwB(fcL, fcR, norder, f), ss = s;
			for(int i=1; i<norderT; i++) ss *= s;
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

	void Smooth(const float* sig, float* sigout, const int n, int half_l, bool abs = true) {
		/* copy and compute abs of the signal into sigw */
		int i, j, wb, we;
		long double wsum;
		const float *sigw = sig;
		std::unique_ptr<float> sigw_p;
		if( abs ) {
			sigw_p.reset(new float[n]);
			float *sigwp = sigw_p.get(); 
			for( i=0; i<n; i++ ) sigwp[i] = fabs(sig[i]);
			sigw = sigwp;
		}

		//
		if(half_l*2>n-1) half_l = (n-1)/2;
		for(i=0,wsum=0.;i<=half_l;i++) wsum += sigw[i];
		wb = 0; we = i;
		for(i=1;i<=half_l;i++,we++) {
			sigout[i-1] = wsum/we;
			wsum += sigw[we];
		}
	   for(j=we;i<n-half_l;i++,wb++,we++) {
			//if( i>80000/dt && i<82000/dt ) std::cerr<<(i-1)*dt<<" "<<sig[i-1]<<" "<<wsum<<" / "<<j<<std::endl;
			sigout[i-1] = wsum/j;
	      wsum += ( sigw[we] - sigw[wb] );
		}
	   for(;i<n;i++,wb++) {
			sigout[i-1] = wsum/(we-wb);
			wsum -= sigw[wb];
	   }
		if(wsum>1.e-15) {
			sigout[n-1] = wsum/(we-wb);
		}
	}

	void FillGap( float *pbeg, float *pend, float mean1, float mean2, float std, int hlen, int step ) {
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
			//for(int i=1; i<step; i++) *(p+i) = *p + slope * i;
			for(int i=1; i<step; i++) *(p-step+i) = *(p-step) + slope * i;  // by lili
      }
		// generate tapered random numbers for the 2nd half
		for(; p<pend; p+=step) { 
			int ndiff = (pend-p);
			float gdamp = exp( alpha * ndiff * ndiff );
			*p = rand() * gdamp; 
			slope = ( *p - *(p-step) ) * oostep;
			//for(int i=1; i<step; i++) *(p+i) = *p + slope * i;
         for(int i=1; i<step; i++) *(p-step+i) = *(p-step) + slope * i;	// by lili
		}
		// connect the last several points
		p -= step; if( p == pend-1 ) return;
		auto plast = pend-1; *plast = rand();
		slope = ( *plast - *p ) / ( plast - p );
      for(p=p+1; p<plast; p++) *p = *(p-1) + slope;
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
      //float ftmp = 2./(ns*dt);
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

	int read_rec(int rec_flag, const char *fname, int len, int *rec_b, int *rec_e, int *nrec) {
		FILE *frec;
		int irec;
		if( rec_flag ) {
			if((frec = fopen(fname,"r")) == NULL) return 0;
			//pthread_mutex_lock(&fiolock);
			#pragma omp critical(sacIO)
			for(irec=0;;irec++)
				if(fscanf(frec,"%d %d", &rec_b[irec], &rec_e[irec])!=2) break;
			*nrec=irec;
			fclose(frec);
			//pthread_mutex_unlock(&fiolock);
			if(irec==0) return 0;
		}
		else {
			rec_b[0]=0; rec_e[0]=len-1;
			*nrec=1;
		}
		return 1;
	}
	void UpdateRec(const std::string& name, int *rec_b, int *rec_e, int nrec) {
		int rec_b1[1000], rec_e1[1000], nrec1;
		FILE *frec;
		int irec, irec1;
		if( ! read_rec(1, name.c_str(), 0, rec_b1, rec_e1, &nrec1) ) {
			std::cerr<<"*** Warning("<<FuncName<<"): cannot open record file "<<name<<" ***"<<std::endl;
			frec = fopen(name.c_str(), "w");
			for(irec=0; irec<nrec; irec++)
				fprintf(frec, "%d %d\n", rec_b[irec], rec_e[irec]);
			fclose(frec);
		}
		int recB, recE;
		std::string name2 = name + "2";
		frec = fopen(name2.c_str(), "w");
		for(irec=0; irec<nrec; irec++)
			for(irec1=0;irec1<nrec1;irec1++){
				if(rec_b[irec]>=rec_e1[irec1]) continue;
				if(rec_e[irec]<=rec_b1[irec1]) break;
				recB = std::max(rec_b[irec],rec_b1[irec1]);
				recE = std::min(rec_e[irec],rec_e1[irec1]);
				fprintf(frec, "%d %d\n", recB, recE);
			}
		fclose(frec);
	}

	#ifdef DISAZI_H
	void ComputeDisAzi( SAC_HD& shd ) {
		try {
			Path<float> path( shd.evlo, shd.evla, shd.stlo, shd.stla );
			shd.dist = path.Dist(); shd.az = path.Azi1(); 
			shd.baz = path.Azi2(); shd.baz += shd.baz<180.?180.:-180.;
		} catch (const std::exception& e) {}
	}
	#else
	void ComputeDisAzi( SAC_HD& shd ) { throw ErrorSR::UndefMethod( FuncName, "DisAzi.h not included!" ); }
	#endif

	// stockwell transform
	//char *Wisfile = NULL;
	//char *Wistemplate = "%s/.fftwis";
	//#define WISLEN 8
	std::string Wisfile;

	void set_wisfile(void) {
		if (! Wisfile.empty()) return;
		char *home;	home = getenv("HOME");
		Wisfile = std::string(home) + "/.fftwis";
		//Wisfile = (char *)malloc(strlen(home) + WISLEN + 1);
		//sprintf(Wisfile, Wistemplate, home);
	}

	/* Convert frequencies in Hz into rows of the ST, given sampling rate and length. */

	int st_freq(double f, int len, double srate)	{
		return floor(f * len / srate + .5);
	}

	/* Stockwell transform of the real array data. The len argument is the
		number of time points, and it need not be a power of two. The lo and hi
		arguments specify the range of frequencies to return, in Hz. If they are
		both zero, they default to lo = 0 and hi = len / 2. The result is
		returned in the complex array result, which must be preallocated, with
		n rows and len columns, where n is hi - lo + 1. For the default values of
		lo and hi, n is len / 2 + 1. */
	int planlen = 0;
	double *g;
	fftw_plan p1, p2;
	fftw_complex *h, *H, *G;
	void st(int len, int lo, int hi, float *data, double *result)	{
		int i, k, n, l2;
		double s, *p;
		FILE *wisdom;

		/* Check for frequency defaults. */

		if (lo == 0 && hi == 0) {
			hi = len / 2;
		}

		/* Keep the arrays and plans around from last time, since this
			is a very common case. Reallocate them if they change. */

		if (len != planlen && planlen > 0) {
			#pragma omp critical(fftw)
	      {
			fftw_destroy_plan(p1);
			fftw_destroy_plan(p2);
			}
			fftw_free(h);
			fftw_free(H);
			fftw_free(G);
			free(g);
			planlen = 0;
		}

		if (planlen == 0) {
			planlen = len;
			h = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * len);
			H = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * len);
			G = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * len);
			g = (double *)malloc(sizeof(double) * len);

			/* Get any accumulated wisdom. */

			set_wisfile();
			#pragma omp critical (fftw)
			{
			wisdom = fopen(Wisfile.c_str(), "r");
			if (wisdom) {
				fftw_import_wisdom_from_file(wisdom);
				fclose(wisdom);
			}

			/* Set up the fftw plans. */
			p1 = fftw_plan_dft_1d(len, h, H, FFTW_FORWARD, FFTW_MEASURE);
			p2 = fftw_plan_dft_1d(len, G, h, FFTW_BACKWARD, FFTW_MEASURE);

			/* Save the wisdom. */
			wisdom = fopen(Wisfile.c_str(), "w");
			if (wisdom) {
				fftw_export_wisdom_to_file(wisdom);
				fclose(wisdom);
			}
			}
		}

		/* Convert the input to complex. Also compute the mean. */

		s = 0.;
		memset(h, 0, sizeof(fftw_complex) * len);
		for (i = 0; i < len; i++) {
			h[i][0] = data[i];
			s += data[i];
		}
		s /= len;

		/* FFT. */

		fftw_execute(p1); /* h -> H */

		/* Hilbert transform. The upper half-circle gets multiplied by
			two, and the lower half-circle gets set to zero.  The real axis
			is left alone. */

		l2 = (len + 1) / 2;
		for (i = 1; i < l2; i++) {
			H[i][0] *= 2.;
			H[i][1] *= 2.;
		}
		l2 = len / 2 + 1;
		for (i = l2; i < len; i++) {
			H[i][0] = 0.;
			H[i][1] = 0.;
		}

		/* Fill in rows of the result. */

		p = result;

		/* The row for lo == 0 contains the mean. */

		n = lo;
		if (n == 0) {
			for (i = 0; i < len; i++) {
				*p++ = s;
				*p++ = 0.;
			}
			n++;
		}

		/* Subsequent rows contain the inverse FFT of the spectrum
			multiplied with the FFT of scaled gaussians. */

		while (n <= hi) {

			/* Scale the FFT of the gaussian. Negative frequencies
				wrap around. */

			g[0] = gauss(n, 0);
			l2 = len / 2 + 1;
			for (i = 1; i < l2; i++) {
				g[i] = g[len - i] = gauss(n, i);
			}

			for (i = 0; i < len; i++) {
				s = g[i];
				k = n + i;
				if (k >= len) k -= len;
				G[i][0] = H[k][0] * s;
				G[i][1] = H[k][1] * s;
			}

			/* Inverse FFT the result to get the next row. */

			fftw_execute(p2); /* G -> h */
			for (i = 0; i < len; i++) {
				*p++ = h[i][0] / len;
				*p++ = h[i][1] / len;
			}

			/* Go to the next row. */

			n++;
		}
	}

	/* This is the Fourier Transform of a Gaussian. */

	static double gauss(int n, int m)
	{
		return exp(-2. * M_PI * M_PI * m * m / (n * n));
	}

	/* Inverse Stockwell transform. */

	void ist(int len, int lo, int hi, const double *data, float *result, 
				const int nskipb = 0, const int nskipe = 0)	{
		int i, n, l2;
		FILE *wisdom;

		/* Check for frequency defaults. */

		if (lo == 0 && hi == 0) {
			hi = len / 2;
		}

		/* Keep the arrays and plans around from last time, since this
			is a very common case. Reallocate them if they change. */

		if (len != planlen && planlen > 0) {
			#pragma omp critical(fftw)
	      {
			fftw_destroy_plan(p2);
			}
			fftw_free(h);
			fftw_free(H);
			planlen = 0;
		}

		if (planlen == 0) {
			planlen = len;
			h = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * len);
			H = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * len);

			/* Get any accumulated wisdom. */

			set_wisfile();
			#pragma omp critical (fftw)
			{
			wisdom = fopen(Wisfile.c_str(), "r");
			if (wisdom) {
				fftw_import_wisdom_from_file(wisdom);
				fclose(wisdom);
			}

			/* Set up the fftw plans. */
			p2 = fftw_plan_dft_1d(len, H, h, FFTW_BACKWARD, FFTW_MEASURE);

			/* Save the wisdom. */

			wisdom = fopen(Wisfile.c_str(), "w");
			if (wisdom) {
				fftw_export_wisdom_to_file(wisdom);
				fclose(wisdom);
			}
			}
		}

		/* Sum the complex array across time. */

		memset(H, 0, sizeof(fftw_complex) * len);
		const double *p = data;
		for (n = lo; n <= hi; n++) {
			for (i = 0; i < len; i++) {
				H[n][0] += *p++;
				H[n][1] += *p++;
			}
		}

		/* Invert the Hilbert transform. */

		l2 = (len + 1) / 2;
		for (i = 1; i < l2; i++) {
			H[i][0] /= 2.;
			H[i][1] /= 2.;
		}
		l2 = len / 2 + 1;
		for (i = l2; i < len; i++) {
			H[i][0] = H[len - i][0];
			H[i][1] = -H[len - i][1];
		}

		/* Inverse FFT. */

		fftw_execute(p2); /* H -> h */
		float *pres = result+nskipb;
		for (i = nskipb; i < len-nskipe; i++) {
			*pres++ = h[i][0] / len;
		}
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
	shd(recin.shd), pimpl( new SRimpl(*(recin.pimpl)) ) {
	if( recin.sig && shd.npts>0 ) {
		sig.reset(new float[shd.npts]);
		if( ! sig )
			throw ErrorSR::MemError( FuncName, "new failed!");
		std::copy(recin.sig.get(), recin.sig.get()+recin.shd.npts, sig.get()); 
	}
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
	if( recin.sig && npts>0 ) {
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
   sig.reset(new float[recin.shd.npts]() );
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
   std::ifstream fsac(fname);
   if( ! fsac )
		throw ErrorSR::BadFile( FuncName, "reading from " + fname );
   //if( SHDMap.empty() ) pimpl->CreateSHDMap();
   //pthread_mutex_lock(&fiolock);
	size_t rdsize = sizeof(SAC_HD);
	#pragma omp critical(sacIO)
   fsac.read( reinterpret_cast<char *>(&shd), rdsize );
	if( fsac.gcount() != rdsize )
		throw ErrorSR::BadFile( FuncName, "failed to retrieve sac header from "+fname+
													 "("+std::to_string(fsac.gcount())+std::to_string(rdsize)+")" );
   fsac.close();
   //pthread_mutex_unlock(&fiolock);
	// make sure shd.e is filled
	if( shd.npts > 0 ) shd.e = shd.b + shd.delta * (shd.npts-1);
}

/* load sac header+signal from file 'fname', memory is allocated on heap */
void SacRec::Load () {
   std::ifstream fsac(fname);
   if( ! fsac )
		throw ErrorSR::BadFile( FuncName, "reading from " + fname );
   //if( SHDMap.empty() ) pimpl->CreateSHDMap();
   //pthread_mutex_lock(&fiolock);
	size_t rdsize = sizeof(SAC_HD);
	#pragma omp critical(sacIO)
   fsac.read( reinterpret_cast<char *>(&shd), rdsize );
	if( fsac.gcount() != rdsize )
		throw ErrorSR::BadFile( FuncName, "failed to retrieve sac header from " + fname );
   //pthread_mutex_unlock(&fiolock);
	// make sure shd.e is filled
	if( shd.npts > 0 ) shd.e = shd.b + shd.delta * (shd.npts-1);
	// allocate memory for sac signal
	if( shd.npts <= 0 )
		throw ErrorSR::BadParam( FuncName, "negative npts("+std::to_string(shd.npts)+") in header");
   sig.reset(new float[shd.npts]);
	if( ! sig )
		throw ErrorSR::MemError( FuncName, "new failed!");
   float* sigsac = sig.get();
	rdsize = sizeof(float) * shd.npts;
	#pragma omp critical(sacIO)
   fsac.read( reinterpret_cast<char *>(sigsac), rdsize );
	if( fsac.gcount() != rdsize )
		throw ErrorSR::BadFile( FuncName, "failed to extract sac sig from " + fname );
   fsac.close();

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

void SacRec::Dump( const std::string fname, float tb, float te ) const {
	if( !sig || shd.npts<=0 )
		throw ErrorSR::EmptySig(FuncName);
	bool tofile = !fname.empty();
	std::ofstream fout(fname);
	if( tofile && !fout )
		throw ErrorSR::BadFile( FuncName, "writing to " + fname );
	std::ostream& sout = tofile ? fout : std::cout;
	int ib = tb==NaN ? 0 : Index(tb);
	int ie = te==NaN ? shd.npts : Index(te);
	float b = shd.b, dt = shd.delta, *sigsac = sig.get();
	for( int i=ib; i<ie; i++ )
		sout<<b+i*dt<<" "<<sigsac[i]<<"\n";
}

/* dump signal to a vector of PointCs */
void SacRec::Dump( std::vector<PointC>& dataV, float tb, float te ) const {
	if( !sig || shd.npts<=0 )
		throw ErrorSR::EmptySig(FuncName);
	float b = shd.b, dt = shd.delta, *sigsac = sig.get();
	int ib = tb==NaN ? 0 : Index(tb);
	int ie = te==NaN ? shd.npts : Index(te);
	dataV.clear(); dataV.reserve(ie-ib);
	for( int i=ib; i<ie; i++ )
		dataV.push_back( PointC(b+i*dt, sigsac[i]) );
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
void SacRec::sqrt() {
   if( !sig )
		throw ErrorSR::EmptySig(FuncName);
	float *sigsac = sig.get();
   for(int i=0; i<shd.npts; i++) sigsac[i] = std::sqrt(sigsac[i]);
}

void SacRec::Mul( const float mul ) {
   if( !sig )
		throw ErrorSR::EmptySig(FuncName);
	float *sigsac = sig.get();
   for(int i=0; i<shd.npts; i++) sigsac[i] *= mul;
}

void SacRec::Add( const float val ) {
   if( !sig ) throw ErrorSR::EmptySig(FuncName);
	float *sigsac = sig.get();
   for(int i=0; i<shd.npts; i++) 
		sigsac[i] += val;
}

void SacRec::Addf( const SacRec& sac2 ) {
   if( ! sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( ! sig ) { *this = sac2; return; }
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );
	//const auto& sig2 = sac2.sig;
	float *sigsac = sig.get(), *sigsac2 = sac2.sig.get();
   for(int i=0; i<shd.npts; i++) 
		sigsac[i] += sigsac2[i];
	if(shd.user0!=NaN && sac2.shd.user0!=NaN) shd.user0 += sac2.shd.user0;
}

void SacRec::Subf( const SacRec& sac2 ) {
   if( !sig || !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( ! sig ) { *this = sac2; this->Mul(-1.); return; }
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );
	//const auto& sig2 = sac2.sig;
	float *sigsac = sig.get(), *sigsac2 = sac2.sig.get();
   for(int i=0; i<shd.npts; i++) 
		sigsac[i] -= sigsac2[i];
}

void SacRec::Mulf( const SacRec& sac2 ) {
   if( !sig || !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );
	//const auto& sig2 = sac2.sig;
	float *sigsac = sig.get(), *sigsac2 = sac2.sig.get();
   for(int i=0; i<shd.npts; i++) 
		if( sigsac2[i]!=0. ) sigsac[i] *= sigsac2[i];
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
float SacRec::BAzi() const {
	auto shdl = shd;
	if( shdl.az == NaN ) pimpl->ComputeDisAzi( shdl );
	return shdl.baz;
}
float SacRec::BAzi() {
	if( shd.az == NaN ) pimpl->ComputeDisAzi( shd );
	return shd.baz;
}

double SacRec::DayTime() const {
	return 3600.*shd.nzhour + 60.*shd.nzmin + shd.nzsec + 0.001*shd.nzmsec;
}

int SacRec::AbsDay(int year0) const {
	if( year0 > shd.nzyear )
		throw ErrorSR::BadParam(FuncName, "year0("+std::to_string(year0)+") > shd.nzyear("+std::to_string(shd.nzyear)+")");
   int nyday = 0;
   for( int i=year0+1; i<shd.nzyear; i++ ) {
      if ( (i%400==0) || (i%100!=0&&i%4==0) ) nyday += 366;
      else nyday += 365;
   }
	return nyday + shd.nzjday;
}

double SacRec::AbsTime(int year0) const {
   //if( ! sig ) return -1.;
   //if( shd == sac_null ) return -1.; // operator== not defined yet
   if( shd.npts <= 0 ) return -1.;
   if( shd.nzjday == NaN || shd.nzyear == NaN || shd.nzhour == NaN ||
       shd.nzmin == NaN || shd.nzsec == NaN || shd.nzmsec == NaN ) return -1;
   //computes time in s relative to year0 (defaulted to 1900)
   return 24.*3600.*AbsDay() + DayTime();
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
   rms = std::sqrt(rms/(neff-1));
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

float SacRec::PrecNoiseRatio() const {
	float noise = 0., cor_pre = 0.;
	int i, ndis, nb;
	int lag = (shd.npts-1)/2;
	float dt=shd.delta;
	float dist=shd.dist;
	ndis = (int)floor((dist/0.8+50.)/dt+0.5);
	if( ndis > lag ) return NaN;
	nb = int(lag*4/5);
	if( ndis > nb ) nb = ndis;
	auto sigsac = sig.get();
	for( i = lag+nb; i< 2*lag; i++) noise += sigsac[i]*sigsac[i];
	noise = std::sqrt(noise/(lag-nb));

	ndis = (int)floor((dist/4.5-50.)/dt+0.5);
	if( ndis < 10./dt ) return NaN;
	nb = int(100./dt);
	if(ndis < nb ) nb = ndis;
	for( i = lag-nb; i<= lag+nb; i++) cor_pre += sigsac[i]*sigsac[i];
	cor_pre = std::sqrt(cor_pre/(2.*nb));
	return cor_pre/noise;
}

/* compute accurate sig value at a given time (fit with a parabola) */
float SacRec::Sig( double time ) const {
	int itime = Index(time);
#ifndef PARABOLA_H
	return sig.at(itime);
#else
	if( itime == 0 ) itime += 1;
	else if( itime == shd.npts-1 ) itime -= 1;
	float* sigsac = sig.get();
	// to avoid rounding problem, shift to zero time
	PointC p1(-shd.delta, sigsac[itime-1]);
	PointC p2(0.,   sigsac[itime]  );
	PointC p3(shd.delta, sigsac[itime+1]);
	//std::cerr<<p1<<"\n"<<p2<<"\n"<<p3<<"\nresult: "<<time<<" "<<Parabola( p1, p2, p3 )(time)<<std::endl;
	return Parabola( p1, p2, p3 )(time-X(itime));
#endif
}

/* smoothing (running average) */
void SacRec::Smooth( float timehlen, SacRec& sacout, bool abs, float tb, float te ) const {
   if( ! sig )
		throw ErrorSR::EmptySig(FuncName);

   int half_l = (int)floor(timehlen/shd.delta+0.5);
	if( half_l < 1 ) {
		sacout = *this; return;
	}

	//std::cerr<<"SacRec::Smooth: new npts = "<<shd.npts<<std::endl;
	if( tb!=NaN || te!=NaN ) { // copy sac into sacout if partialy smoothing{
		sacout = *this;
	} else {	// resize sacout.sig otherwise
		sacout.shd = shd;
		sacout.sig.reset( new float[shd.npts] );
		if( ! sacout.sig )
			throw ErrorSR::MemError( FuncName, "new failed!");
	}

	float *sigsac = sig.get();
   float *sigout = sacout.sig.get();

	int npts = te==NaN ? shd.npts : Index(te);
	if( tb != NaN ) {
		int nb = Index(tb); npts -= nb;
		sigsac += nb; sigout += nb;
	}
	//std::cerr<<"SacRec::Smooth: before call "<<tb<<" "<<npts<<" "<<half_l<<"   "<<sigsac<<std::endl;
	pimpl->Smooth(sigsac, sigout, npts, half_l, abs);
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
   std = std::sqrt(std/(dataV.size()-1));

	return true;
}

float SacRec::MeanPha(const float fb, const float fe) const {
	if( ! sig )
		throw ErrorSR::EmptySig(FuncName);
	int ib = fb==NaN ? 0 : Index(fb);
	int ie = fe==NaN ? shd.npts : Index(fe);
	float *sigph = sig.get();
	double x = 0., y = 0.;
	for(int i=ib; i<ie; i++) {
		x += cos(sigph[i]);
		y += sin(sigph[i]);
	}
	return atan2(y, x);
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
void SacRec::ToAmPh( SacRec& sac_am, SacRec& sac_ph, const float fl, const float fu, const int nfout ) const {
	if( shd.npts > maxnpts4parallel ) {
		#pragma omp critical(largesig)
		{
		ToAmPh_p(sac_am, sac_ph, fl, fu, nfout);
		}
	} else {
		ToAmPh_p(sac_am, sac_ph, fl, fu, nfout);
	}
}
void SacRec::ToAmPh_p( SacRec& sac_am, SacRec& sac_ph, const float fl, const float fu, const int nfout ) const {
   if( ! sig ) 	// check signal
		throw ErrorSR::EmptySig(FuncName, fname);
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
	float deltaf = 1./(delta*ns);
	if( amp==nullptr || pha==nullptr )
		throw ErrorSR::MemError( FuncName, "new failed!");
	if( fl>=0. && fu>0. && fu>fl ) {
		for(int i=0; i<nk; i++) {
			fftw_complex& cur = sf[i];
			float Abtw = pimpl->btwB(fl, fu, 8, deltaf*i);
			amp[i] = std::sqrt(cur[0]*cur[0] + cur[1]*cur[1]) * delta * Abtw;
			pha[i] = atan2(cur[1], cur[0]);
		}
	} else {
		for(int i=0; i<nk; i++) {
			fftw_complex& cur = sf[i];
			amp[i] = std::sqrt(cur[0]*cur[0] + cur[1]*cur[1]) * delta;
			pha[i] = atan2(cur[1], cur[0]);
		}
	}
   //fftw_free(s); 
	fftw_free(sf);
   sac_am.sig.reset(amp);
	sac_ph.sig.reset(pha);

   sac_am.shd.npts = nk;
   sac_am.shd.delta = deltaf;
   sac_am.shd.b = 0.;
	sac_am.shd.e = sac_am.shd.delta * (nk-1);
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
	sac_re.shd.e = sac_re.shd.delta * (nk-1);
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
void SacRec::Filter ( double f1, double f2, double f3, double f4, const int type, SacRec& srout, bool zeroPhase ) const {
	if( shd.npts > maxnpts4parallel ) {
		#pragma omp critical(largesig)
		{
		Filter_p(f1, f2, f3, f4, type, srout, zeroPhase);
		}
	} else {
		Filter_p(f1, f2, f3, f4, type, srout, zeroPhase);
	}
}
void SacRec::Filter_p ( double f1, double f2, double f3, double f4, const int type, SacRec& srout, bool zeroPhase ) const {
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

/* Stockwell Transform and its inverse */
std::vector<double> SacRec::SWT( int& ifl, int& ifu, int& itb, int& ite,
											float fl, float fu, float tb, float te ) const {
	if( !sig || shd.npts<=0 )
		throw ErrorSR::EmptySig(FuncName);
	// correct time&freq boundaries
	const float deltat = shd.delta;
	const float fmax = 0.5 / deltat;
	if( fl < 0. ) fl = 0.;
   if( fu<=fl || fu>fmax ) fu = fmax;
	const float tmax = shd.b + deltat*(shd.npts-1);
	if( tb<shd.b || tb==NaN ) tb = shd.b;
	if( te<=tb || te>tmax || te==NaN ) te = tmax;
	// indexes
	itb = Index(tb); ite = Index(te);
	int nt = ite - itb;
	const float deltaf = 1 / (deltat*nt); 
	// extend freq range for the btw taper
	ifl = nint(fl / deltaf); ifu = nint(fu / deltaf);
	const float toler = 0.01;
/*
	int btworder = 6;
	for(; ifl>0; ifl--) {
		float freq = ifl * deltaf;
		if( pimpl->btwB(fl, fu, btworder, freq) < toler ) break;
	}
	for(; ifu<nint(fmax/deltaf); ifu++) {
		float freq = ifu * deltaf;
		if( pimpl->btwB(fl, fu, btworder, freq) < toler ) break;
	}
*/
	// allocate output vector
	int nf = ifu-ifl+1;
	std::vector<double> datastV;
	try {
		datastV.resize(2*nt*nf);
	} catch (const std::exception& e) {
		throw ErrorSR::MemError( FuncName, "SWT vector alloc failed!");
	}
	// call st
	pimpl->st(nt, ifl, ifu, &sig[itb], &datastV[0]);
	// apply taper to datastV
/*
	auto Idatast = datastV.begin();
	for(int ifreq=ifl; ifreq<=ifu; ifreq++) { // on freq
		float freq = ifreq * deltaf;
		float btw = pimpl->btwB(fl, fu, btworder, freq);
		if( btw != btw ) btw = 0.;
		for(int itime=0; itime<nt; itime++) { // on time
			double &real = *(Idatast++);
			double &imag = *(Idatast++);
			real *= btw; imag *= btw;
		}
	}
*/
	return datastV;
}

// zero out any signal < 
void SacRec::NoiseZeroOut( SacRec& sacout, std::vector<int>& recb, std::vector<int>& rece, const float tlen_min,
									const float nofactor, const float nomin, const float ttaper ) const {
	if( &sacout == this )
		throw ErrorSR::BadParam(FuncName, "ZeroOut in place is not allowed!");
	// smooth the original signal
	Smooth(tlen_min, sacout);
	// and find the reference noise level
	auto sbeg = &(sacout.sig[0]);
	int inoiselevel = (int)(shd.npts*0.5);
	std::nth_element(sbeg, sbeg+inoiselevel, sbeg+shd.npts);
	float noiselevel = sacout.sig[inoiselevel];
	// signal < (noiselevel*nfactor or nomin) is not allowed
	float stdmin = std::max(nomin, noiselevel * nofactor);
	// identify invalid (near zero) segments as defined by stdmin
	sacout = *this;
	int len_min = nint(tlen_min/shd.delta), npsame = 1;
	const float *sigsac = sig.get();
	float *sigosac = sacout.sig.get();
	recb.clear(); rece.clear();
	recb.push_back(0);
	float sigmin = sigsac[0], sigmax = sigsac[0];
	int npts = sacout.shd.npts, nrece = npts;
	for(int i=1; i<npts; i++) {
		auto &sigcur = sigsac[i];
		if( sigcur > sigmax ) sigmax = sigcur;
		else if( sigcur < sigmin ) sigmin = sigcur;
		if( sigmax-sigmin>stdmin || i==npts-1 ) {
			if( npsame >= len_min ) {
				int ib = i-npsame, ie = i==npts-1 ? npts : i;
				float tsameb = sacout.X(ib);
				sacout.cosTaperR(tsameb-ttaper, tsameb, false);
				float tsamee = sacout.X(ie-1);
				sacout.cosTaperL(tsamee, tsamee+ttaper, false);
				for(int j=i-npsame; j<ie; j++) sigosac[j] = 0.;
				if( ib == 0 ) {
					recb[0] = ie;
				} else if( ie == npts ) {
					nrece = ib;
				} else {
					rece.push_back(ib);
					recb.push_back(ie);
				}
			}
			npsame = 1;
			sigmin = sigmax = sigcur;
		} else {
			npsame++;
		}
	}
	rece.push_back(nrece);
	if( recb[0] == 0 ) sacout.cosTaperL(shd.b, shd.b+ttaper);
	if( nrece == npts ) sacout.cosTaperR(shd.e-ttaper, shd.e);
}


void SacRec::ISWT( const std::vector<double>& datastV, 
						 const int ifl, const int ifu, const int itb, const int ite,
						 const int nskipb, const int nskipe ) {
	if( !sig || shd.npts<=0 )
		throw ErrorSR::EmptySig(FuncName);
	int len = ite - itb;
	const float deltat = shd.delta, deltaf = 1 / (deltat*len);
	// check if datastV matches the given parameters
	int n = ifu-ifl+1;
	if( datastV.size() != 2*len*n )
		throw ErrorSR::BadParam(FuncName, "datastV size does not match the input parameters");
	// call ist
	pimpl->ist(len, ifl, ifu, &datastV[0], &sig[itb], nskipb, nskipe);
}

/* ---------------------- t-f normalization with stockwell transform ---------------------- */
// this is a 2-D one bit normalization in the f-t domain, stablized with the cutoff_factor:
// cutoff_factor = 0 (cut any point>the_smallest) - 1 (cut any point>the_largest --not modified at all)
void SacRec::STNorm( SacRec& sacout, const float thlen, float fl, float fu, 
							float tsafe, const float t_seg, const short norm_method ) const {
	if( !sig || shd.npts<=0 )
		throw ErrorSR::EmptySig(FuncName);

	if( &sacout == this )
		throw ErrorSR::BadParam(FuncName, "STNorm in place is not allowed!");
	sacout.MutateAs(*this);

	if( tsafe < thlen )
		std::cerr<<"Warning("<<FuncName<<"): tsafe smaller than thlen"<<std::endl;

	// correct freq boundaries
	const float deltat = shd.delta;
	const float fmax = 0.5 / deltat;
	if( fl < 0. ) fl = 0.;
   if( fu<=fl || fu>fmax ) fu = fmax;

	// zero out near-zero time windows to stablize the normalization
	SacRec sac_0(*this);
	NoiseZeroOut( sac_0 );

	// divide the signal into small segments and perform stockwell transform 
	// on each of them for better time/memory performance
	float ttotal = (shd.npts-1) * deltat;
	float timespan = ttotal - tsafe*2.;
	int nseg = std::max(1, nint(timespan / t_seg));
	float tseg = timespan / nseg;
	int nskip = (int)floor(tsafe/deltat);

	// loop over segments
	const float cutoff_refperc  = 0.6, cut_factor = 3.;
	for(int iseg=0; iseg<nseg; iseg++) {
		// stockwell transform
		// extend each ends by tsafe to prevent edge effect
		// so that the extended segment becomes (safe _ data _ safe)
		float tbsafe = shd.b + iseg*tseg;
		float tb = tbsafe + tsafe;
		float te = tb + tseg + shd.delta;
		float tesafe = te + tsafe;
		int ifl, ifu, itb, ite;
		std::vector<double> datastV = sac_0.SWT( ifl, ifu, itb, ite, fl, fu, tbsafe, tesafe );
		int nf = ifu-ifl+1, nt = ite-itb;
		const float deltaf = 1 / (deltat*nt); 

/*
std::ofstream fout("debug1.txt");
auto Idatast1 = datastV.begin();
for(int ifreq=ifl; ifreq<ifu; ifreq++) { // on freq
	float freq = ifreq * deltaf;
	for(int itime=0; itime<nt; itime++) { // on time
		float time = (itb+itime) * deltat + shd.b;
		double &real = *(Idatast1++);
		double &imag = *(Idatast1++);
		double amp = sqrt(real*real+imag*imag);
		fout<<freq<<" "<<time<<" "<<amp<<"\n";
	}
	fout<<"\n\n";
}
fout.close(); fout.clear();
*/
		// normalize at each freq
		auto sigsac0 = sac_0.sig.get();
		auto Idatast = datastV.begin();
		for(int ifreq=ifl; ifreq<=ifu; ifreq++) { // on freq
			float freq = ifreq * deltaf;
			auto Ib = Idatast;
			// compute/save amps of the current freq slice
			std::vector<float> dataV_fslice(nt);
			auto *dfslice = &dataV_fslice[0];
			for(int itime=0; itime<nt; itime++) { // on time
				double &real = *(Idatast++);
				double &imag = *(Idatast++);
				double amp = std::sqrt(real*real+imag*imag);
				dfslice[itime] = amp;
			}
			// smooth the current time slide
			int half_l = (int)floor(thlen/shd.delta+0.5);
			std::vector<float> dataV_fslice_s(nt);
			auto *dfslice_s = &dataV_fslice_s[0];
			pimpl->Smooth(dfslice, dfslice_s, nt, half_l);
			// and perform selected normalization
			Idatast = Ib;
			if( norm_method == 0 ) {
				for(int itime=0; itime<nt; itime++, Idatast+=2) {
					//if( dfslice[itime] == 0 ) continue;
					double &real = *(Idatast);
					double &imag = *(Idatast+1);
					float mul = 1. / dfslice_s[itime];
					if( sigsac0[itb+itime] == 0. ) mul = 1.0e-10;
					real *= mul; imag *= mul;
				}
			} else if( norm_method == 1 ) {
				//// find the nth smallest point where n=cutoff_factor*nt
				// for each freq, find the cutoff_refval at the nth smallest point where n=cutoff_refperc*npts
				// and press down anything > cutoff_refval * cut_factor
				int ip_cut = floor(cutoff_refperc * nt);
				if( ip_cut >= nt ) ip_cut = nt - 1;
				std::vector<float> dataV_tmp(dataV_fslice);
				std::nth_element(dataV_tmp.begin(), dataV_tmp.begin()+ip_cut, dataV_tmp.end());
				float amp_ref = dataV_tmp[ip_cut];
				float amp_refS = amp_ref * amp_ref;
				//float crct = 1., f2end = (ifu-ifreq) * deltaf;
				//if( f2end < 0.03 ) crct = pow(f2end/0.03, 100.);
				// loop again, cut any points > amp_cut
				for(int itime=0; itime<nt; itime++, Idatast+=2) {
					double &real = *(Idatast);
					double &imag = *(Idatast+1);
					float mul;
					if( sigsac0[itb+itime] == 0. ) {
						mul = 1.0e-10;	// set zero if sig_o == 0
					} else {
						//float amp = dfslice[itime];
						float amps = dfslice_s[itime];
						//mul = amp_cut / amps;
						mul = (12*amps*amps+amp_refS) / ((8*amps*amps+5*amp_refS)*amps);
					}
					real *= mul; imag *= mul;
				}
			} else {
				throw ErrorSR::BadParam(FuncName, "unknown norm method");
			}
		}

/*
fout.open("debug2.txt");
Idatast1 = datastV.begin();
for(int ifreq=ifl; ifreq<ifu; ifreq++) { // on freq
	float freq = ifreq * deltaf;
	for(int itime=0; itime<nt; itime++) { // on time
		float time = (itb+itime) * deltat + shd.b;
		double &real = *(Idatast1++);
		double &imag = *(Idatast1++);
		double amp = sqrt(real*real+imag*imag);
		fout<<freq<<" "<<time<<" "<<amp<<"\n";
	}
	fout<<"\n\n";
}
fout.close(); fout.clear();
*/
		// inverse transform
		// only write from tb to te and discard anything outside
		int nskipb = iseg==0 ? 0 : nskip;
		int nskipe = iseg==nseg-1 ? 0 : nskip;
//std::cerr<<tbsafe<<" "<<tb<<" "<<te<<" "<<tesafe<<"   "<<nskipb<<" "<<nskipe<<std::endl;
		sacout.ISWT( datastV, ifl, ifu, itb, ite, nskipb, nskipe );
		auto sigout = sacout.sig.get();
		for(int i=itb; i<ite; i++) 
			if( sigsac0[i] == 0. ) sigout[i] = 0.;
	}

	/*
	// to reduce edge effects from SWT
	// (1) taper at the two ends
	float tend = shd.b + ttotal;
	float ttaper = tsafe * 1.5;
	sacout.cosTaperL( shd.b, shd.b+ttaper );
	sacout.cosTaperR( tend-ttaper,   tend );
	// (2) and around any zero segments
	// smooth the original signal
	float tlen_min = 30.;
	SacRec sac_sm; Smooth(tlen_min, sac_sm);
	// and find the reference noise level
	auto sbeg = &(sac_sm.sig[0]);
	int inoiselevel = (int)(shd.npts*0.5);
	std::nth_element(sbeg, sbeg+inoiselevel, sbeg+shd.npts);
	float noiselevel = sac_sm.sig[inoiselevel];
	float stdmin = noiselevel * 0.1;
	// identify invalid (near zero) segments as defined by stdmin
	int len_min = nint(tlen_min/shd.delta), npsame = 1;
	const float *sigsac = sig.get();
	float *sigosac = sacout.sig.get();
	float siglast = sigsac[0];
	for(int i=1; i<sacout.shd.npts; i++) {
		if( fabs(sigsac[i]-siglast) < stdmin ) {
			npsame++;
		} else {
			if( npsame >= len_min ) {
				float tsameb = sacout.X(i-npsame);
				sacout.cosTaperR(tsameb-ttaper, tsameb, false);
				float tsamee = sacout.X(i-1);
				sacout.cosTaperL(tsamee, tsamee+ttaper, false);
				for(int j=i-npsame; j<i; j++) sigosac[j] = 0.;
//std::cerr<<tsameb<<" - "<<tsamee<<" taperred!"<<std::endl;
			}
			npsame = 1;
		}
		siglast = sigsac[i];
	}
	*/
}


/* cosine tapper */
void SacRec::cosTaperL( const float fl, const float fh, bool zeroout ) {
	if( !sig || shd.npts<=0 )	// check signal
		throw ErrorSR::EmptySig(FuncName);

	const float& delta = shd.delta;
	int i, il = (int)ceil((fl-shd.b)/delta);
	float *sigsac = sig.get();
	if( zeroout ) for( i=0; i<il; i++ ) sigsac[i] = 0.;
	else i = std::max(il, 0);
	float finit = shd.b+i*delta, fwidth = fh-fl, dfh = fh-finit;
	float ftmp = M_PI / fwidth;
   for(; dfh>0&&i<shd.npts; i++, dfh-=delta) {
		float amplif = ( 1. + cos(ftmp*dfh) ) * 0.5;
		sigsac[i] *= amplif;
   }

   return;
}
void SacRec::cosTaperR( const float fl, const float fh, bool zeroout ) {
   if( !sig || shd.npts<=0 )	// check signal
		throw ErrorSR::EmptySig(FuncName);

	const float& delta = shd.delta;
   int i = std::max(0, (int)ceil((fl-shd.b)/delta));
	float finit = shd.b+i*delta, fwidth = fh-fl, dfl = finit-fl;
	float ftmp = M_PI / fwidth;
	float *sigsac = sig.get();
   for(; i<shd.npts&&dfl<fwidth; i++, dfl+=delta) {
		float amplif = ( 1. + cos(ftmp*dfl) ) * 0.5;
		sigsac[i] *= amplif;
   }
   if( zeroout ) for(;i<shd.npts;i++) sigsac[i] = 0.;

   return;
}

void SacRec::btwTaperB( double fcL, double fcR, int norder ) {
   if( !sig || shd.npts<=0 )	// check signal
		throw ErrorSR::EmptySig(FuncName);
	auto sigsac = sig.get();
	double f = 0.; const auto &delta = shd.delta;
	for( int i=0; i<shd.npts; i++, f+=delta )
		sigsac[i] *= pimpl->btwB(fcL, fcR, norder, f);
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
void SacRec::cut( float tb, float te, SacRec& sac_result ) const {
	if( !sig || shd.npts<=0 )	// check signal
		throw ErrorSR::EmptySig(FuncName);

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


bool SacRec::merge( SacRec sacrec2 ) {
   // make sure that both signals are loaded
   if( !sig && !sacrec2.sig )
		throw ErrorSR::EmptySig(FuncName, "both sacs are empty");

	if( !sig || !sacrec2.sig ) {
		if( ! sig ) *this = sacrec2;
		return false;
	}
   
   SAC_HD& shd2 = sacrec2.shd;
   // starting and ending time
   double t1b = this->AbsTime() + shd.b;
   double t2b = sacrec2.AbsTime() + sacrec2.shd.b;
   double t1e = t1b + (shd.npts-1)*shd.delta;
   double t2e = t2b + (shd2.npts-1)*shd2.delta;

	// merge in place if sacrec2 is within sacrec1
   double dt = (int)floor(shd.delta*1e8+0.5)/1e8;
   int nb = (int)floor((t2b-t1b)/dt+0.5);
	double tshift = (nb*dt-(t2b-t1b));// + (shd2.b-shd.b);
   if( fabs(tshift) > 1.e-3 )
		(*report) << "signal shifted by " << tshift << "secs" << std::endl;
	if( t2b>=t1b && t2e<=t1e ) {
		std::copy(&(sacrec2.sig[0]), &(sacrec2.sig[0])+shd2.npts, &(sig[nb]));
		return true;
	}

	// merge out of place otherwise
   double T1 = std::min(t1b, t2b), T2 = std::max(t1e, t2e);

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
   bool reversed = false;
   if( t1b > t2b ) {
      reversed = true;
      nb = (int)floor((t1b-t2b)/dt+0.5);
		tshift = (nb*dt-(t1b-t2b));// + (shd.b-shd2.b);
   }

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

	return true;
}

int SacRec::arrange(const char* recname) {
   // count holes
   float maxfloat = std::numeric_limits<float>::max()*0.99;
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

void SacRec::ZoomToEvent( const SAC_HD& eshd, float evlon, float evlat, float tb, float tlen, std::string ename ) {
	if( evlon == NaN ) evlon = eshd.evlo;
	if( evlat == NaN ) evlat = eshd.evla;
	if( tb == NaN ) tb = eshd.b;
	if( tlen == NaN ) tlen = (eshd.npts-1) * eshd.delta;
	if( ename.empty() ) ename = std::to_string(eshd.nzyear) + "-" + std::to_string(eshd.nzjday) + "-" +
										 std::to_string(eshd.nzhour) + "-" + std::to_string(eshd.nzmin) + "-" + std::to_string(eshd.nzsec);
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
	pimpl->ComputeDisAzi( shd );
}

/* remove mean and trend */
void SacRec::RTrend() {
   // fit a*x+b with at most ~10000 points
   int npts = shd.npts;
	int step = std::max(1, nint(npts/10000.));
	// fit starts
	float *sigsac = sig.get();
   float X = 0., Y = 0., X2 = 0., Y2 = 0., XY = 0.;
	int nptse = 0;
   for(int i=0;i<npts;i+=step) {
      X += i;
      Y += sigsac[i];
      X2 += i*i;
      Y2 += sigsac[i]*sigsac[i];
      XY += i*sigsac[i];
		nptse++;
   }
   float a = (nptse*XY-X*Y)/(nptse*X2-X*X);
   float b = (-X*XY+X2*Y)/(nptse*X2-X*X);
	if( a!=a || b!=b )
		throw ErrorSR::FloatOverflow(FuncName, "fitting params");
   // correct sig and DEPMEN
   float mean = 0., max = sigsac[0], min = sigsac[0];
   float shift = b;
   for(int i=0;i<npts;i++,shift+=a) {
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
	int hourmid = (int)ceil( (DayTime() + shd.b) / 3600. ); 
   sprintf(buff, "%s %s %s %4d %3d %f %f %d -t %d: -f %s -v >& /dev/null", evrexe.c_str(), sta, ch, shd.nzyear, shd.nzjday, f1, f4, nf, hourmid, fresp.c_str());
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
	//std::cerr<<fresp<<" "<<nameam<<" "<<nameph<<std::endl;
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
		throw ErrorSR::ExternalError( FuncName, std::to_string(sizeA) + " AMP file(s) found from " + fresp);
	if( ! openA )
		throw ErrorSR::BadFile( FuncName, "reading from " + std::string(nameam) );
	if( sizeP != 1 )
		throw ErrorSR::ExternalError( FuncName, std::to_string(sizeP) + " PHASE file(s) found from " + fresp);
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
//Write("debug2.SAC");
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
	if( nptst < 1 ) throw ErrorSR::InsufData(FuncName, "npts<=0 after resampling");
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

void SacRec::Interpolate( int npts_ratio, SacRec& sac2 ) const {
	if( npts_ratio <= 1 ) return;
	if( ! sig )
		throw ErrorSR::EmptySig(FuncName);
	
	// prepare sac2
	int npts2 = (shd.npts-1) * npts_ratio + 1;
	sac2.shd = shd;
	sac2.shd.npts = npts2; sac2.shd.delta /= npts_ratio;
	sac2.sig.reset();
	sac2.sig.reset( new float[npts2] );

	// interpolate starts
	auto sigsac = sig.get(), sig2sac = sac2.sig.get();
	for(int i1=0; i1<shd.npts-1; i1++) {
		float sig0 = sigsac[i1];
		float siginc = (sigsac[i1+1] - sigsac[i1]) / npts_ratio;
		int i2 = i1*npts_ratio;
		sig2sac[i2] = sig0;
		for(int i=1; i<npts_ratio; i++) 
			sig2sac[i2+i] = sig0 + i*siginc;
	}
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

void RunAvg( float timehlen, float Eperl, float Eperh, std::vector<SacRec>& sacV, bool normByFirst ) {
	SacRec sac_sigmax;
	// denominator (maximum smoothed signal)
	for( auto &sac : sacV ) {
		SacRec sac_sm;
		/* filter into the earthquake band */
		if( Eperl != -1. ) {
			SacRec sac_eqk;
			float f2 = 1./Eperh, f1 = f2*0.6, f3 = 1./Eperl, f4 = f3*1.4;
			sac.BandpassCOSFilt( f1, f2, f3, f4, sac_eqk );
			sac_eqk.Smooth( timehlen, sac_sm );
		} else {
			sac.Smooth( timehlen, sac_sm );
		}
		sac_sigmax.PullUpTo( sac_sm );
		if( normByFirst ) break;
	}
	// apply normalizer
	for( auto &sac : sacV ) sac.Divf(sac_sigmax);
}

void SacRec::Whiten( float fl, float fu, float fhlen ) {
	SacRec sac_am, sac_ph;
	ToAmPh(sac_am, sac_ph);
	sac_am.RunAvg(fhlen, -1, -1);
	sac_am.btwTaperB(fl, fu);
	FromAmPh(sac_am, sac_ph);
}

// earthquake cutting

bool SacRec::EqkCut( SacRec& sacout, std::vector<int>& rec_b, std::vector<int>& rec_e, 
							const float Eperl, const float Eperu, const float maxnoise_factor, 
							bool apptaper, const std::string& recname ) const {
	// any point with val<=sigmin is assumed 0
	float sigmin = 1.;
	// time length of taper window to reduce FFT edge effects
   float ttaper = 200.;
	int nptaper = nint(ttaper / shd.delta);
	// evenly sample 1000 points. assume invalid if >60% are zeros
	auto sigsac = sig.get();
   int n = shd.npts, ninc = n>1000?n/1000:1, npole=0;
   for(int i=0;i<n;i+=ninc) if(fabs(sigsac[i])<sigmin) npole++;
   if(npole>600) {
      std::cerr<<"*** Warning("<<FuncName<<"): Signal time length not long enough. ***"<<std::endl;
      return false;
   }
	// apply eqk filter for event detection
   //float* sigw = new float[n];
   //double f2 = 1./Eperh, f1 = f2*0.8, f3 = 1./Eperl, f4 = f3*1.2;
	SacRec sacw;
	int ib = 0, ie = n;
   if( Eperl == -1 ) {
		sacw = *this;
	} else {
		BandpassBTWFilt( 1./Eperu, 1./Eperl, 6, sacw );
		// apply taper on both side to reduce edge effects
		for(ib=0; ib<shd.npts&&sigsac[ib]<sigmin; ib++){} float tb = X(ib);
		sacw.cosTaperL( tb, tb+ttaper );
		for(ie=shd.npts-1; ie>0&&sigsac[ie]<sigmin; ie--){} float te = X(ie);
		sacw.cosTaperR( te-ttaper, te );
	}
	float* sigw = sacw.sig.get();

   // noise window npts (twindow sec length)
   double dt = (double)(shd.delta), twindow = 1000.;
   int s1k=(int)floor(twindow/dt+0.5); // npts of a twindow sec window
   int nos1k = (int)(n/s1k);

   // for (each of) the ith twindow sec window, search for maximum amplitude and store into win_max[i]
   double win_max[nos1k];
   memset (win_max,0,nos1k*sizeof(double));
   for(int i=0;i<n;i++) sigw[i] = fabs(sigw[i]);
   int ii, is;
   for( ii=0,is=0; is<nos1k; is++ ){
      for( ; ii<(is+1)*s1k; ii++ ) 
         if(win_max[is]<sigw[ii]) win_max[is]=sigw[ii]; 
   }
   for( ;ii<n;ii++ )
      if(win_max[is]<sigw[ii]) win_max[is]=sigw[ii];

   // sort win_max
   std::vector<double> win_max_sorted( win_max, win_max+nos1k );
   std::sort( win_max_sorted.begin(), win_max_sorted.end() );

	// discard any window with max<=sigmin (which, with an unit of nm, is pretty much 0)
   std::vector<double>::iterator itermin;
   for(itermin=win_max_sorted.begin(); itermin<win_max_sorted.end(); itermin++) if( *itermin > sigmin ) break;
	/*
   // and define max noise level as 3 x average_of_the_smallest_20_windows
   double noisemax = 0.;
   if( itermin < win_max_sorted.end() ) {
      for(auto iter=itermin; iter<win_max_sorted.end() && iter<itermin+20; iter++) noisemax += *iter; 
      noisemax *= 3./(iter-itermin);
   }
	*/
	// and define max noise level as 3 x median (assuming no more than 40% of the record is affected by eqks)
	int iwinb = itermin-win_max_sorted.begin(), imid = (nos1k+iwinb)*0.6;
	double noisemax = win_max_sorted[imid] * 3.;
	//std::cerr<<"noise-max = "<<noisemax<<std::endl;

   // compute noise average and noise std between windows
   double window_avg = 0.;
	auto iter = itermin;
	while( iter<win_max_sorted.end() && *iter<noisemax ) window_avg += *(iter++);
   ii = iter-itermin;
   window_avg /= ii;
   double window_std=0., dtmp;
   for(iter=itermin; iter<itermin+ii; iter++) {
      dtmp = window_avg-*iter;
      window_std += dtmp * dtmp;
   }
   window_std=std::sqrt(window_std/(ii-1));
	//std::cerr<<window_avg<<" "<<window_std<<std::endl;

   // mark windows with a max amp > window_avg+maxnoise_factor*window_std to be 'zero'
   dtmp = window_avg+maxnoise_factor*window_std;
   short keep[nos1k];
   for( int i =0; i<nos1k; i++) keep[i] = win_max[i] > dtmp ? 0 : 1;

   // and zero out invalidated windows
	sacout = *this;
	auto sigosac = sacout.sig.get();
   for( int i=0; i < nos1k; i++)
      if( keep[i] == 0 ) {
			//std::cerr<<"zero out points "<<i*s1k<<" to "<<(i+1)*s1k<<std::endl;
			for( ii=i*s1k; ii<(i+1)*s1k; ii++) sigosac[ii] = 0.;
		}

   // locate contigious valid windows and apply cosine tapers
	float ttaper_sig = apptaper ? ttaper : 0.;
	sacout.NoiseZeroOut(rec_b, rec_e, 30., 0.1, sigmin, ttaper_sig);

	// remove invalid windows with a length < 2.5*twindow sec
 	int winlen_min = (int)ceil(2.5*twindow/dt);
	sigosac = sacout.sig.get();	// re-point! NoiseZeroOut invalidated the old sigosac pointer
	for(int irec=0; irec<rec_b.size(); ) {
		if( rec_e[irec]-rec_b[irec] > winlen_min ) {	irec++; continue;	}
		for(int i=rec_b[irec]; i<rec_e[irec]; i++) sigosac[i] = 0.;
		rec_b.erase(rec_b.begin()+irec);
		rec_e.erase(rec_e.begin()+irec);
	}
	//for(int irec=0; irec<rec_b.size(); irec++) std::cerr<<"   "<<X(rec_b[irec])<<" "<<X(rec_e[irec])<<std::endl;

	/*
   // locate all rec_begin and rec_end pairs, zero out the windows that are shorter than winlen_min
	int rec_i=0;
	rec_b.reserve(50); rec_e.reserve(50);
	rec_b.resize(1); rec_e.resize(1);
   rec_b[0] = 0;
   for( int i=1; i<nos1k; i++){ 
      if(keep[i]-keep[i-1] == 1) rec_b[rec_i]=i*s1k; // a new window begins
      else if(keep[i]-keep[i-1] == -1) { // the current window ends
         rec_e[rec_i]=i*s1k;
			// invalidate the current window if it is shorter than winlen_min
         if ((rec_e[rec_i]-rec_b[rec_i]) < winlen_min) {
            for(ii=rec_b[rec_i]; ii<rec_e[rec_i]; ii++) sig[ii] = 0.;
         } else {
				rec_i++; 
				rec_b.resize(rec_i+1); rec_e.resize(rec_i+1);
			}
      }
   }
   // mark the last rec_end and check its window length 
   if(keep[nos1k-1]==1) {
      if( (n-rec_b[rec_i]) < winlen_min*0.6 )
         for(ii=rec_b[rec_i]; ii<n; ii++) sig[ii]=0;
      else { rec_e[rec_i] = n; rec_i++; }
   }
	*/

   // check if there's enough data (20%) left
	ii = 0;
   for(int i=0; i<rec_b.size(); i++) ii += rec_e[i] - rec_b[i];
   if( ii < 0.2*n ) {
      std::cerr<<"*** Warning("<<FuncName<<"): Time length <20\% after removing earthquakes ("<<fname<<") ***";
      return false;
   }

	/*
   // taper ttaper sec of data on each end of the windows to be safe 
   for(int i=0;i<rec_i;i++) {
      if(rec_b[i]!=0){
			float tb = X(rec_b[i]-1);
			cosTaperL( tb, tb+ttaper, false );
		}
		rec_b[i] += 0.5*nptaper;
      if(rec_e[i]!=n){
			float te = X(rec_e[i]+1);
			cosTaperR( te-ttaper, te, false );
		}
      rec_e[i] -= 0.5*nptaper;
   }
	*/

   // produce a new rec file named ft_name_rec2
   if( ! recname.empty() ) pimpl->UpdateRec(recname, &rec_b[0], &rec_e[0], rec_b.size());
   /* norm by running average if required */
   //if( tnorm_flag == 4 ) RunAvgNorm( sig, shd, sigw );

   return true;
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

void CalcTransferF( const SacRec& sac1, const SacRec& sac2, const float fmin, SacRec& Coh, SacRec& Adm, SacRec& Pha ) {
	SacRec sac1_am, sac1_ph; sac1.ToAmPh(sac1_am, sac1_ph);
	SacRec sac2_am, sac2_ph; sac2.ToAmPh(sac2_am, sac2_ph);
	//sac1_am.Mul(1.0e-5);	sac2_am.Mul(1.0e-5); // to prevent floating-point overflow
	// calculate autospectral density functions Gss, Grr and
	// the one-sided cross-spectral density function Grs
	SacRec Grr(sac1_am);  Grr.Mulf(sac1_am);
	SacRec Gss(sac2_am);  Gss.Mulf(sac2_am);
	SacRec GrsA(sac2_am); GrsA.Mulf(sac1_am);
	SacRec GrsP(sac2_ph); GrsP.Subf(sac1_ph);
	// convert GrsAmp&GrsPha to GrsR&GrsI in place
   AmPhToReIm( GrsA, GrsP );
	// freqency domain smoothing for statistical stability
	int nsm = 20; float fhlen = nsm * GrsA.shd.delta;
	//std::cerr<<fhlen<<std::endl;
	//Grr.Write("Grr.SAC"); Gss.Write("Gss.SAC"); GrsA.Write("GrsA.SAC"); GrsP.Write("GrsP.SAC");
	GrsA.Smooth(fhlen, false, fmin); GrsP.Smooth(fhlen, false, fmin);
	Grr.Smooth(fhlen, true, fmin); Gss.Smooth(fhlen, true, fmin);
	//Grr.Write("Grr.SAC_sm"); Gss.Write("Gss.SAC_sm"); GrsA.Write("GrsA.SAC_sm"); GrsP.Write("GrsP.SAC_sm");
	// convert GrsR&GrsI back to GrsAmp&GrsPha
	ReImToAmPh( GrsA, GrsP );
	// construct coherence, admittance, and phase from the spectral density functions
	Pha = std::move(GrsP); Adm = GrsA; Adm.Divf( Gss );
	Coh = Adm; Coh.Mulf( GrsA ); Coh.Divf( Grr ); Coh.sqrt();
	//Coh.Write("debugCoh.SAC"); Adm.Write("debugAdm.SAC"); Pha.Write("debugPha.SAC"); 
	// smooth Coh
	//Coh.Smooth(0.002, false, fmin); 
}

void SacRec::CCFromAmPh(SacRec& sac_am, SacRec& sac_ph, const SAC_HD& shd1, const SAC_HD& shd2) {
	SacRec sac_CC;
	sac_CC.FromAmPh(sac_am, sac_ph);
	// form final result
	int lag = shd1.npts-1, ns = (sac_am.shd.npts-1) * 2;
	shd = shd1;
	shd.npts = lag*2 + 1;
	float Tshift = shd2.b-shd1.b, Tlag = lag * shd.delta;
	shd.b = -Tlag + Tshift;
	shd.e = Tlag + Tshift;
	shd.user0 = 1;
	strncpy(shd.kevnm, shd1.kstnm, 8); strncpy(shd.kstnm, shd2.kstnm, 8);
	shd.evlo = shd1.stlo; shd.stlo = shd2.stlo;
	shd.evla = shd1.stla; shd.stla = shd2.stla;
	pimpl->ComputeDisAzi( shd );
	sig.reset( new float[lag*2+1]() );
	float *cor = sig.get(), *CCout = sac_CC.sig.get();
   for( int i = 1; i< (lag+1); i++) {
      cor[lag+i] =  CCout[i];
      cor[lag-i] =  CCout[ns-i];
   }
   cor[lag] = CCout[0];

	// normalize
	//float *outsig = sacout.sig.get();
   //for(int i=0; i<sacout.npts; i++) outsig[i] *= ;
}

SacRec SacRec::CrossCorrelate( SacRec& sac2, const std::string& outname, int ctype ) {
	// check input
	if( !sig || !sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	if( shd.delta != sac2.shd.delta )
		throw ErrorSR::HeaderMismatch(FuncName, "dt");
	if( shd.npts != sac2.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(shd.npts)+" - "+std::to_string(sac2.shd.npts) );

	// FFT on sig1
	SacRec sac1_am, sac1_ph;
	ToAmPh(sac1_am, sac1_ph);
	
	// FFT on sig2
	SacRec sac2_am, sac2_ph;
	sac2.ToAmPh(sac2_am, sac2_ph);

	SacRec saco_am, saco_ph;
	CrossCorrelateSACs( sac1_am, sac1_ph, sac2_am, sac2_ph, saco_am, saco_ph, ctype );

	if(! outname.empty()) {
		saco_am.Write(outname+"_am");
		saco_ph.Wrap(); saco_ph.Write(outname+"_ph");
	}

	// whiten
	//saco_am.RunAvg(0.001, -1, -1);
	//saco_am.cosTaperL(0.01, 0.0125);
	//saco_am.cosTaperR(0.4, 0.5);
	// convert back to time domain (with length ns)
	SacRec sacout; sacout.CCFromAmPh( saco_am, saco_ph, shd, sac2.shd );
	return sacout;
}

SacRec CrossCorrelateSACs( const SacRec& sac1_am, const SacRec& sac1_ph, const SacRec& sac2_am, const SacRec& sac2_ph, 
									const SAC_HD& shd1, const SAC_HD& shd2, int ctype ) {
	SacRec saco_am, saco_ph;
	CrossCorrelateSACs( sac1_am, sac1_ph, sac2_am, sac2_ph, saco_am, saco_ph, ctype );
	SacRec sacout; sacout.CCFromAmPh( saco_am, saco_ph, shd1, shd2 );
	return sacout;
}

/* Cross-Correlate with another sac record
	ctype=0: Cross-Correlate (default) 
	ctype=1: deconvolve (sac.am/sac2.am)
	ctype=2: deconvolve (sac2.am/sac.am) */
void CrossCorrelateSACs( const SacRec& sac1_am, const SacRec& sac1_ph, const SacRec& sac2_am, const SacRec& sac2_ph, 
							SacRec& saco_am, SacRec& saco_ph, int ctype ) {
	// check input
	auto &shd = sac1_am.shd;
	if( ! (sac1_am.sig&&sac1_ph.sig && sac2_am.sig&&sac2_ph.sig) )
		throw ErrorSR::EmptySig(FuncName);
	if( shd.delta != sac2_am.shd.delta )
		throw ErrorSR::HeaderMismatch(FuncName, "dt");
	if( shd.npts != sac2_am.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(shd.npts)+" - "+std::to_string(sac2_am.shd.npts) );

	// initialize amp and phase out
	int npts = shd.npts;
	if( ctype == 2 ) saco_am = sac2_am;
	else saco_am = sac1_am;
	/*
	saco_am.sig.reset( new float[ns]() );
	if( ! saco_am.sig )
		throw ErrorSR::MemError( FuncName, "new failed for saco_am!");
	saco_am.shd = sac1_am.shd;
	saco_am.shd.npts = ns;
	*/

	// correlate in freq domain
	float *amout = saco_am.sig.get();
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
	//sac1_am.clear(); sac2_am.clear();

	// compute phase out
	saco_ph = sac2_ph;
	/*
	saco_ph.sig.reset( new float[ns]() );
	if( ! saco_ph.sig )
		throw ErrorSR::MemError( FuncName, "new failed for saco_ph!");
	saco_ph.shd = sac1_ph.shd;
	saco_ph.shd.npts = ns;
	*/
	float *phout = saco_ph.sig.get();
	float *ph1 = sac1_ph.sig.get(), *ph2 = sac2_ph.sig.get();
   for(int i=0; i<npts; i++) phout[i] -= ph1[i];
	saco_ph.Wrap();

	// clear sacs & free memories 2
	//sac1_ph.clear(); sac2_ph.clear();

}


void AmPhToReIm( SacRec& sac_am, SacRec& sac_ph ) {
	if( ! sac_am.sig&&sac_ph.sig )
		throw ErrorSR::EmptySig(FuncName);
	const int npts = sac_am.shd.npts;
	if( npts != sac_ph.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(npts)+" - "+std::to_string(sac_ph.shd.npts) );
	auto sigam = sac_am.sig.get(), sigph = sac_ph.sig.get();
	for(int i=0; i<npts; i++) {
		float &am = sigam[i], &ph = sigph[i];
		float re = am * cos(ph);
		float im = am * sin(ph);
		am = re; ph = im;
	}
}

void ReImToAmPh( SacRec& sac_re, SacRec& sac_im ) {
	if( ! sac_re.sig&&sac_im.sig )
		throw ErrorSR::EmptySig(FuncName);
	const int npts = sac_re.shd.npts;
	if( npts != sac_im.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(npts)+" - "+std::to_string(sac_re.shd.npts) );
	auto sigre = sac_re.sig.get(), sigim = sac_im.sig.get();
	for(int i=0; i<npts; i++) {
		float &re = sigre[i], &im = sigim[i];
		float am = std::sqrt(im*im + re*re);
		float ph = atan2(im, re);
		re = am; im = ph;
	}
}

void SACRotate( SacRec& sac1, SacRec& sac2, const float deg ) {
	if( fabs(deg) < 1.e-5 ) return;
	if( ! sac1.sig&&sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	const int npts = sac1.shd.npts;
	if( npts != sac2.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(npts)+" - "+std::to_string(sac2.shd.npts) );

	const float degrad = deg * M_PI / 180.;
	const float cosdeg = cos(degrad), sindeg = sin(degrad);
	auto sigsac1 = sac1.sig.get(), sigsac2 = sac2.sig.get();
	for(int i=0; i<npts; i++) {
		float sig1 = sigsac1[i], sig2 = sigsac2[i];
		sigsac1[i] = sig1 * cosdeg + sig2 * sindeg;
		sigsac2[i] = sig2 * cosdeg - sig1 * sindeg;
	}
}

SacRec SACProject( const SacRec& sac1, const SacRec& sac2, const float deg ) {
	if( fabs(deg) < 1.e-5 ) return sac1;
	if( ! sac1.sig&&sac2.sig )
		throw ErrorSR::EmptySig(FuncName);
	const int npts = sac1.shd.npts;
	if( npts != sac2.shd.npts )
		throw ErrorSR::HeaderMismatch(FuncName, "npts: "+std::to_string(npts)+" - "+std::to_string(sac2.shd.npts) );

	const float degrad = deg * M_PI / 180.;
	const float cosdeg = cos(degrad), sindeg = sin(degrad);
	SacRec sacout; sacout.MutateAs(sac1);
	auto sigsac1 = sac1.sig.get(), sigsac2 = sac2.sig.get();
	sacout.Transform2i( [&](const int i, float& val){
		val = sigsac1[i] * cosdeg + sigsac2[i] * sindeg;
	} );
	return sacout;
}


void DumpSACs( const std::vector<SacRec>& sacV, const std::string& outname ) {
	if( sacV.empty() ) return;
	// get header info
	auto &shd0 = sacV[0].shd;
	auto npts = shd0.npts;
	auto b = shd0.b, delta = shd0.delta;
	// get pointers to each sac signal (to help performance)
	float *sigsacA[sacV.size()]; int nvalid = 1;
	sigsacA[0] = sacV[0].sig.get();
	// and check header consistency at the same time
	for(int i=1; i<sacV.size(); i++) {
		if( ! sacV[i].sig ) continue;
		sigsacA[nvalid++] = sacV[i].sig.get();
		const auto &shdc = sacV[i].shd;
		if( shdc.npts!=npts || shdc.delta!=delta || shdc.b!=b )
			throw ErrorSR::HeaderMismatch(FuncName, sacV[i].fname+" npts||delta||b");
	}
	// open out file
	std::ofstream fout(outname);
	if( ! fout ) throw std::runtime_error("IO failed on file "+outname);
	// output
	for( int i=0; i<npts; i++ ) {
		float x = b + delta * i;
		fout<<x<<"\t";
		for( int isac=0; isac<nvalid; isac++ ) fout<<sigsacA[isac][i]<<" ";
		fout<<"\n";
	}
}

