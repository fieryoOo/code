#ifndef SACREC_H
#define SACREC_H

#include "mysac64.h"
#include "MyOMP.h"
#include "Parabola.h"
//#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <limits>
#include <stdexcept>
#include <cmath>

#ifndef FuncName
#define FuncName __FUNCTION__
#endif
/*
// a workaround for nullptr when compiled by old gcc compiler
const class {										// this is a const object...
public:
   template<class T>								// convertible to any type
   operator T*() const { return 0 }			// of null non-member pointer...
   template<class C, class T>					// or any type of null
   operator T C::*() const { return 0; }	// member pointer...
private:
  void operator&() const;						// whose address can't be taken
} nullptr = {};									// and whose name is nullptr
*/


namespace ErrorSR {

   class Base : public std::runtime_error {
   public:
		Base( const std::string funcname, const std::string message )
			: runtime_error(funcname + ": " + message) {
				//PrintStacktrace();
			}
	};

   class BadFile : public Base {
   public:
      BadFile(const std::string funcname, const std::string info = "")
         : Base(funcname, "Invalid or non-accessable file ("+info+").") {}
   };

   class BadParam : public Base {
   public:
      BadParam(const std::string funcname, const std::string info = "")
         : Base(funcname, "Bad parameters ("+info+").") {}
   };

   class EmptySig : public Base {
   public:
      EmptySig(const std::string funcname, const std::string info = "")
         : Base(funcname, "No sac signal loaded in the memory ("+info+").") {}
   };

   class InsufData : public Base {
   public:
      InsufData(const std::string funcname, const std::string info = "")
         : Base(funcname, "Insufficient data points ("+info+").") {}
   };

   class HeaderMismatch : public Base {
   public:
      HeaderMismatch(const std::string funcname, const std::string info = "")
         : Base(funcname, "Header mismatch ("+info+").") {}
   };

   class FloatOverflow : public Base {
   public:
      FloatOverflow(const std::string funcname, const std::string info = "")
         : Base(funcname, "Floating-point overflow ("+info+").") {}
   };

   class UndefMethod : public Base {
   public:
		UndefMethod(const std::string funcname, const std::string info = "")
         : Base(funcname, "Undefined method ("+info+").") {}
   };

   class ExternalError : public Base {
   public:
		ExternalError(const std::string funcname, const std::string info = "")
         : Base(funcname, "External error ("+info+").") {}
   };

   class MemError : public Base {
   public:
		MemError(const std::string funcname, const std::string info = "")
         : Base(funcname, "Memory error ("+info+").") {}
   };
};


//enum SACTYPE { TIME, AMP, PHA, DISP };

class SacRec {
public:
   std::string fname;			// input file name
   SAC_HD shd;				// sac header
   std::unique_ptr<float[]> sig;	// pointer to the signal
   //std::auto_ptr<float> sig;
public:
   /* ------------------------------ con/destructors and operators ------------------------------ */
   /* constructors */
	SacRec( std::ostream& reportin = std::cerr );	// default 1
   SacRec( const std::string& fnamein, std::ostream& reportin = std::cerr );	// default 2
	SacRec( const size_t npts, std::ostream& reportin = std::cerr );	// default 3
   SacRec( const SacRec& recin );		// copy
   SacRec( SacRec&& recin );			// move
   /* operators */
   SacRec &operator= ( const SacRec& recin );	// assignment
   SacRec &operator= ( SacRec&& recin );	// move
   /* destructor */
   ~SacRec(); 

	/* check signal/header validation */
	void updateDeps();
	bool isValid() {
		updateDeps();
		return (shd.delta>0 && shd.npts>0 && shd.depmin==shd.depmin && shd.depmax==shd.depmax);
	}
	bool isZero() {
		float* sigsac = sig.get();
		for(int i=0; i<shd.npts; i++)	if( sigsac[i] != 0. ) return false;
		return true;
	}

	/* assign header and allocate memory */
	void MutateAs ( const SacRec& recin );

	/* allocate memory for signal */
	void ResizeSig() { ResizeSig(shd.npts); }
	void ResizeSig( const size_t npts ) {
		if( npts <= 0 )
			throw ErrorSR::BadParam( FuncName, "negative npts!");
		shd.npts = npts; shd.e = shd.b + shd.delta*npts;
		sig.reset(new float[npts]);
		if( ! sig )
			throw ErrorSR::MemError( FuncName, "new failed!");
	}

   /* ------------------------------ sac file read/write ------------------------------ */
   /* load sac header from file 'fname' */
   void LoadHD();
   void LoadHD( const std::string& fnamein ) { fname = fnamein; LoadHD(); }
   /* read sac header+signal from file 'fname', memory is allocated on heap */
   void Load();
   void Load( const std::string& fnamein ) { fname = fnamein; Load(); }
	/* clear sac and release memory */
	void clear() { sig.reset(); shd = sac_null; fname.clear(); }
   /* write to file '*fname' */
   void WriteHD( const std::string& fname );
   void WriteHD() { WriteHD(fname); }
   void Write( const std::string& fname );
   void Write() { Write(fname); };
	/* load a txt file */
	void LoadTXT( const std::string& fname );
	/* dump signal to stdout/txt */
	void Dump( const std::string fname = "", float tb = NaN, float te = NaN ) const;
	/* dump signal to a vector of PointCs */
	void Dump( std::vector<PointC>& dataV, float tb = NaN, float te = NaN ) const;
	/* dump header to stdout/txt */
	void DumpHD( const std::string fname = "" ) const;
	/* print a single header field */
	void PrintHD( const std::string field, std::ostream &o = std::cout ) const;
   /* change a single header filed */
	void ChHdr(const std::string& field, const std::string& value){
		std::istringstream sin( field + " " + value ); sin >> shd;
	}

	float Dis() const;
	float Dis();
	float Azi() const;
	float Azi();
	float BAzi() const;
	float BAzi();

	const std::string ntname() const { return shd.ntname(); }
	const std::string evname() const { return shd.evname(); }
	const std::string stname() const { return shd.stname(); }
	const std::string chname() const { return shd.chname(); }

   /* ------------------------------ header/signal information ------------------------------ */
	inline size_t Index( const float time, const bool correctOFB = false ) const {
		int index = nint((time-shd.b)/shd.delta);
		if( index<0 || index> shd.npts )
			if( correctOFB ) {
				if( index<0 ) index = 0;
				else index = shd.npts;
			} else {
				throw ErrorSR::BadParam(FuncName, "index out of range");
			}
		return index; 
	}
	inline double X( const size_t index ) const { return (double)shd.b + index*shd.delta; }
   /* compute the absolute time in sec relative to 1900.01.00 */
	int AbsDay(int year0 = 1900) const;
   double AbsTime(int year0 = 1900) const;
	double DayTime() const;
   /* update/reformat header time if shd.nzmsec is modified and is out of the range [0,1000) */
   void UpdateTime();
   /* search for min&max signal positions and amplitudes */
	void MinMax (int& imin, int& imax) const { MinMax( imin, imax, shd.b, shd.e ); }
	void MinMax (int& imin, int& imax, float tbegin, float tend) const;
   void MinMax ( float tbegin, float tend, float& tmin, float& min, float& tmax, float& max ) const;
   /* compute the root-mean-square average in a given window */
   float RMSAvg ( float tbegin, float tend ) const { return RMSAvg( tbegin, tend, 1); }
   float RMSAvg ( float tbegin, float tend, int step ) const;
	bool Mean ( float& mean ) const { return Mean(shd.b, shd.e, mean); }
	bool Mean ( float tbegin, float tend, float& mean ) const { return Mean(tbegin, tend, 1, mean); }
	bool Mean ( float tbegin, float tend, int step, float& mean ) const;
	bool MeanStd ( float& mean, float& std ) const { return MeanStd(shd.b, shd.e, mean, std); }
	bool MeanStd ( float tbegin, float tend, float& mean, float& std ) const { return MeanStd(tbegin, tend, 1, mean, std); }
	bool MeanStd ( float tbegin, float tend, int step, float& mean, float& std ) const;
	float MeanPha( const float fb = NaN, const float fe = NaN ) const;
	// compute accurate time/amplitude of the peak (fit with a parabola)
	float Tpeak() const { float t, a; Peak(t, a); return t; }
	float Tpeak( const float tbegin, const float tend ) const { float t, a; Peak(t, a, tbegin, tend); return t; }
	float Apeak() const { float t, a; Peak(t, a); return a; }
	float Apeak( const float tbegin, const float tend ) const { float t, a; Peak(t, a, tbegin, tend); return a; }
	void Peak(float& tpeak, float& apeak) const { Peak(tpeak, apeak, shd.b, shd.e); }
	void Peak(float& tpeak, float& apeak, const float tbegin, const float tend) const;
	float Sig( double time ) const;	// compute accurate sig value at a given time (fit with a parabola)

	float PrecNoiseRatio() const;

	void NoiseZeroOut( std::vector<int>& recb, std::vector<int>& rece, const float tlen_min = 30., 
							 const float nofactor = 0.1, const float nomin = 5., const float ttaper = 200. ) {
		SacRec sacout; 
		NoiseZeroOut( sacout, recb, rece, tlen_min, nofactor, nomin, ttaper );
		*this = std::move(sacout);
	}
	void NoiseZeroOut( SacRec& sacout, const float tlen_min = 30., 
							 const float nofactor = 0.1, const float nomin = 5., const float ttaper = 200. ) const {
		std::vector<int> recb, rece;
		NoiseZeroOut( sacout, recb, rece, tlen_min, nofactor, nomin, ttaper );
	}
	void NoiseZeroOut( SacRec& sacout, std::vector<int>& recb, std::vector<int>& rece, const float tlen_min = 30., 
							 const float nofactor = 0.1, const float nomin = 5., const float ttaper = 200. ) const;

   /* ------------------------------ single-sac operations ------------------------------ */
	template<class Functor>	void Transform(const Functor& func, size_t ib=0, int ie=NaN) {
		Foreach(func, ib, ie);
	}
	template<class Functor>	void Transform2(const Functor& func, size_t ib=0, int ie=NaN) {
		Foreach2(func, ib, ie);
	}
	template<class Functor>	void Transform2i(const Functor& func, size_t ib=0, int ie=NaN) {
		Foreach2i(func, ib, ie);
	}

	template<class Functor>	void Foreach(const Functor& func, size_t ib=0, int ie=NaN) {
		if( !sig ) throw ErrorSR::EmptySig(FuncName);
		if( ib < 0 ) ib = 0;
		if( ie==NaN || ie>shd.npts ) ie = shd.npts;
		float* sigsac = sig.get();
		for(int i=ib; i<ie; i++)	func(sigsac[i]);
	}
	template<class Functor>	void Foreach2(const Functor& func, size_t ib=0, int ie=NaN) const {
		if( !sig ) throw ErrorSR::EmptySig(FuncName);
		if( ib < 0 ) ib = 0;
		if( ie==NaN || ie>shd.npts ) ie = shd.npts;
		float* sigsac = sig.get();
		for(int i=ib; i<ie; i++)	func(X(i), sigsac[i]);
	}
	template<class Functor>	void Foreach2i(const Functor& func, size_t ib=0, int ie=NaN) const {
		if( !sig ) throw ErrorSR::EmptySig(FuncName);
		if( ib < 0 ) ib = 0;
		if( ie==NaN || ie>shd.npts ) ie = shd.npts;
		float* sigsac = sig.get();
		for(int i=ib; i<ie; i++)	func(i, sigsac[i]);
	}
	template<class Functor>	void Foreach3(const Functor& func, size_t ib=0, int ie=NaN) const {
		if( !sig ) throw ErrorSR::EmptySig(FuncName);
		if( ib < 0 ) ib = 0;
		if( ie==NaN || ie>shd.npts ) ie = shd.npts;
		float* sigsac = sig.get();
		for(int i=ib; i<ie; i++)	func(i, X(i), sigsac[i]);
	}

	void sqrt();
   void Mul( const float mul );
	void Add( const float val );
	void Sub( const float val ) { Add(-val); }
	void Addf( const SacRec& sac2 );
	void Subf( const SacRec& sac2 );
	void Mulf( const SacRec& sac2 );
	void Divf( const SacRec& sac2 );

	void Reverse() { 
		SacRec sac2;
		Reverse(sac2);
		*this = std::move(sac2);
	}
	void Reverse( SacRec& sac2 );

	float SNR( const float tsignall, const float tsignalh, const float tnoisel, const float tnoiseh ) const;

	/* performs integration in the time domain using the trapezoidal rule */
	void IntegrateT() { IntegrateT(*this); }
	void IntegrateT( SacRec& sac_out ) const;

	/* performs integration in the frequency domain (omega arithmetic) */
	void Integrate() { Integrate(*this); }
	void Integrate( SacRec& sac_out ) const;

	/* performs differentiation in the frequency domain (omega arithmetic) */
	void Differentiate() { Differentiate(*this); }
	void Differentiate( SacRec& sac_out ) const;

	/* phase wrap and unwrap */
	void Wrap();
	void Unwrap();

	/* amplitude specturm operations */
	void PullUpTo( const SacRec& sac2 );

	/* FFT & IFFT */
   void ToAm() { ToAm(*this);	}
   void ToAm( SacRec& sac_am ) const {
		SacRec sac_ph;
		ToAmPh( sac_am, sac_ph );
	}
	void FFT( SacRec& sac_re, SacRec& sac_im, const int nfout = 0 ) const;	// in series when sig is large
	void FFT_p( SacRec& sac_re, SacRec& sac_im, const int nfout = 0 ) const;	// always parallel
	void ToAmPh( SacRec& sac_am, SacRec& sac_ph, const float fl=-1., const float fu=-1., const int nfout=0 ) const;	// in series when sig is large
	void ToAmPh_p( SacRec& sac_am, SacRec& sac_ph, const float fl=-1., const float fu=-1., const int nfout=0 ) const;	// always parallel
	void FromAmPh( SacRec& sac_am, SacRec& sac_ph, const short outtype = 0 );		// in series when sig is large
	void FromAmPh_p( SacRec& sac_am, SacRec& sac_ph, const short outtype = 0 );	//	always parallel

	/* Stockwell Transform and its inverse */
	std::vector<double> SWT( float fl = -1., float fu = -1., float tb = NaN, float te = NaN ) const {
		int ifl, ifu, itb, ite; 
		SWT(ifl, ifu, itb, ite, fl, fu, tb, te);
	}
	std::vector<double> SWT( int& ifl, int& ifu, int& itb, int& ite,
									 float fl = -1., float fu = -1., float tb = NaN, float te = NaN ) const;
	void ISWT( const std::vector<double>& datastV, 
				  const int ifl, const int ifu, const int itb, const int te,
				  const int nskipb = 0, const int nskipe = 0 );

	/* filters */
	void LowpassCOSFilt( double f3, double f4 ) { LowpassCOSFilt(f3, f4, *this); }
	void LowpassCOSFilt( double f3, double f4, SacRec& srout ) { Filter(-1., -1., f3, f4, 0, srout); }
	void HighpassCOSFilt( double f1, double f2 ) { HighpassCOSFilt(f1, f2, *this); }
	void HighpassCOSFilt( double f1, double f2, SacRec& srout ) { Filter(f1, f2, -1., -1., 1, srout); }
	void BandpassCOSFilt( double f1, double f2, double f3, double f4 ) { BandpassCOSFilt(f1, f2, f3, f4, *this); }
	void BandpassCOSFilt( double f1, double f2, double f3, double f4, SacRec& srout ) { Filter(f1, f2, f3, f4, 2, srout); }
	void LowpassBTWFilt( double fc, int n ) { LowpassBTWFilt(fc, n, *this); }
	void LowpassBTWFilt( double fc, int n, SacRec& srout ) { Filter(-1., -1., fc, n, 3, srout); }
	void HighpassBTWFilt( double fc, int n ) { HighpassBTWFilt(fc, n, *this); }
	void HighpassBTWFilt( double fc, int n, SacRec& srout ) { Filter(fc, n, -1., -1., 4, srout); }
	void BandpassBTWFilt( double fcL, double fcR, int n, bool zeroPhase = false ) {
		BandpassBTWFilt(fcL, fcR, n, *this, zeroPhase); 
	}
	void BandpassBTWFilt( double fcL, double fcR, int n, SacRec& srout, bool zeroPhase = false ) const { 
		Filter(-1., fcL, fcR, n, 5, srout, zeroPhase); 
		/*
		if( zeroPhase ) {
			srout.Reverse();
			srout.Filter(-1., fcL, fcR, n, 5);
			srout.Reverse();
		} 
		*/
	}
	void GaussianFilt( double fc, double fhlen ) { GaussianFilt(fc, fhlen, *this); }
	void GaussianFilt( double fc, double fhlen, SacRec& srout ) { Filter(-1., fc, fhlen, -1., 6, srout); }
   /* method that performs different types of filters:
	 * type = 0: Lowpass cosine -f3~f4_
	 * type = 1: highpass cosine _f1~f2-
	 * type = 2: bandpass cosine _f1~f2-f3~f4_
	 * type = 3: lowpass butterworth -fc=f3~n=f4_
	 * type = 4: highpass butterworth _fc=f1~n=f2-
	 * type = 5: bandpass butterworth _fcL=f1~nL=f2-fcR=f3~nR=f4_
	 * type = 6: gaussian _fc=f2~fhlen=f3_ */
   void Filter ( double f1, double f2, double f3, double f4, const int type, bool zeroPhase = false ) {	// in-place (in series when sig is large)
		Filter(f1, f2, f3, f4, type, *this, zeroPhase); 
	}
   void Filter ( double f1, double f2, double f3, double f4, const int type, SacRec& srout, bool zeroPhase = false ) const;	// out-of-place (in series when sig is large)
   void Filter_p ( double f1, double f2, double f3, double f4, const int type, SacRec& srout, bool zeroPhase = false ) const;	// always parallel
	/* tapers */
	void cosTaperL( const float fl, const float fh, bool zeroout=true );
	void cosTaperR( const float fl, const float fh, bool zeroout=true );
	void btwTaperB( double fcL, double fcR, int norder = 6 );
	void gauTaper( const float fc, const float fh );
   /* remove mean and trend */
   void RTrend();
   /* remove response and apply filter 
		type: 0=displacement, 1=velocity, 2=acceleration */
   void RmRESP( const std::string& fresp, float perl, float perh, const int type = 1 ) {
		std::string evrexe; RmRESP( fresp, perl, perh, evrexe, type );
	}
   void RmRESP( const std::string& fresp, float perl, float perh, const std::string& evrexe, const int type = 1 );
   /* resample (with anti-aliasing filter) the signal to given sps */
   //void Resample( bool fitParabola = true ) { Resample( floor(1.0/shd.delta+0.5), fitParabola ); }
   void Resample( int sps = -1, bool fitParabola = true );
   void Interpolate( int npts_ratio ) {
		SacRec sac2;
		Interpolate(npts_ratio, sac2);
		*this = std::move(sac2);
	}
   void Interpolate( int npts_ratio, SacRec& sac2 ) const;
	/* smoothing ( running average ) */
	void Smooth( float timehlen, bool abs = true, float tb = NaN, float te = NaN ) { 
		SacRec sac_sm; Smooth(timehlen, sac_sm, abs, tb, te); *this = std::move(sac_sm); 
	}
	void Smooth( float timehlen, SacRec& sacout, bool abs = true, float tb = NaN, float te = NaN ) const;
	void Hilbert() { Hilbert(*this); }
	void Hilbert( SacRec& sacout );
	void Envelope() { Envelope(*this); }
	void Envelope( SacRec& sacout );

   void cut( float tb, float te ) { cut(tb, te, *this); }
   void cut( float tb, float te, SacRec& ) const;
   /* ------------------------------ inter-sac operations ------------------------------ */
   /* merge a second sacrec to the current */
   void Merge( SacRec sacrec2 ) {
      merge( sacrec2 );
      arrange();
   }
   bool merge( SacRec sacrec2 );
   int arrange( const char *recname = nullptr );

	/* ---------- compute the correlation coefficient with an input SacRec ---------- */
	void CCFromAmPh(SacRec& sac_am, SacRec& sac_ph, const SAC_HD& shd1, const SAC_HD& shd2);
	float Correlation( const SacRec& sac2 ) const {	return Correlation( sac2, shd.b, shd.e );	}
	float Correlation( const SacRec& sac2, const float tb, const float te ) const;
	
	friend void CalcTransferF( const SacRec& sac1, const SacRec& sac2, const float fmin, SacRec& Coh, SacRec& Adm, SacRec& Pha );

	/* Cross-Correlate with another sac record
		ctype=0: Cross-Correlate (default) 
		ctype=1: deconvolve (sac.am/sac2.am)
		ctype=2: deconvolve (sac2.am/sac.am) */
	//void CrossCorrelate( SacRec& sac2, const std::string& outname = "" ) { *this = CrossCorrelate(sac2, outname); }
	SacRec CrossCorrelate( SacRec& sac2, const std::string& outname = "", int ctype = 0 );
	friend SacRec CrossCorrelateSACs( const SacRec& sac1_am, const SacRec& sac1_ph, const SacRec& sac2_am, const SacRec& sac2_ph, 
												 const SAC_HD& shd1, const SAC_HD& shd2, int ctype = 0 );
	friend void CrossCorrelateSACs( const SacRec& sac1_am, const SacRec& sac1_ph, const SacRec& sac2_am, const SacRec& sac2_ph, 
											  SacRec& saco_am, SacRec& saco_ph, int ctype = 0 );

	/* ------------------------------- cut by event ---------------------------------- */
	void ZoomToEvent( const std::string etime, float evlon, float evlat, float tb, float tlen, std::string ename = "" );
	void ZoomToEvent( const SAC_HD& eshd, float evlon = NaN, float evlat = NaN, float tb = NaN, float tlen = NaN, std::string ename = "" );

	/* ------------------------------- temporal normalizations ------------------------------- */
	void OneBit();
	void RunAvg( float timehlen, float Eperl, float Eperh );
	friend void RunAvg( float timehlen, float Eperl, float Eperh, std::vector<SacRec>& sacV, bool normByFirst = false );
	void Whiten( float fl, float fu, float fhlen = 0.0002 );
	// EqkCut works on a given seismic record, (attempt to) zero-out all time segments that's 
	// affected by earthquakes or other large amplitude containminations.
	// the percentage of removal is a function of maxnoise_factor
	bool EqkCut( float Eperl=10., float Eperu=40., const float maxnoise_factor = 2.5,
					 bool apptaper=true, const std::string& recname="" ) {
		SacRec sacout;
		std::vector<int> rec_b, rec_e;
		EqkCut( sacout, rec_b, rec_e, Eperl, Eperu, maxnoise_factor, apptaper, recname );
		*this = sacout;
	}
   // mark windows with a max amp > window_avg+maxnoise_factor*window_std to be 'zero'
	bool EqkCut( SacRec& sacout, std::vector<int>& rec_b, std::vector<int>& rec_e, 
					 float Eperl=10., float Eperu=40., const float maxnoise_factor = 2.5, 
					 bool apptaper=true, const std::string& recname="" ) const;

	/* ---------------------- t-f normalization with stockwell transform ---------------------- */
	// this is a 2-D one bit normalization in the f-t domain, stablized with the cutoff_factor:
	// cutoff_factor = 0 (cut any point>the_smallest) - 1 (cut any point>the_largest --not modified at all)
	void STNorm( const float thlen = 40., float fl = -1., float fu = -1., 
					 float tsafe = 200., const float t_seg = 1024., const short nmethod = 0 ) {
		SacRec sacout; STNorm(sacout, thlen, fl, fu, tsafe, t_seg, nmethod);
		*this = sacout;
	}
	void STNorm( SacRec& sacout, const float thlen = 40., float fl = -1., float fu = -1.,
					 float tsafe = 200., const float t_seg = 1024., const short nmethod = 0 ) const;

	/* ------------------------------- memory consumed ------------------------------- */
	float MemConsumed() const;
	void AlwaysParallel() { maxnpts4parallel = std::numeric_limits<int>::max(); }
	// to run the fftw, 16 times the original npts is required ( in&out complex double array with size doubled for specturm ). 20 is used to be safe
	void SetMaxMemForParallel( float MemInMb ) { maxnpts4parallel = (MemInMb * 1024. * 1024. - 1000.) / (4. * 20.); }

   friend void AmPhToReIm( SacRec& sac_am, SacRec& sac_ph );
   friend void ReImToAmPh( SacRec& sac_re, SacRec& sac_im );

   friend void SACRotate( SacRec& sac1, SacRec& sac2, const float deg );
   friend SacRec SACProject( const SacRec& sac1, const SacRec& sac2, const float deg );

   friend void DumpSACs( const std::vector<SacRec>& sacV, const std::string& outname );

	// define string output
   friend std::ostream& operator<< ( std::ostream& o, const SacRec& sac ) {
		o << sac.fname << "  " << sac.evname() << " " << sac.shd.evlo << " " << sac.shd.evla
		  << "  " << sac.stname() << " " << sac.shd.stlo << " " << sac.shd.stla;
		return o;
	}

	static constexpr float NaN = -12345.;

protected:
	int maxnpts4parallel = 1e6;
	static constexpr float twopi = M_PI * 2.0;

private:
   /* impl pointer */
   struct SRimpl;
   std::unique_ptr<SRimpl> pimpl;
	/* reporting stream */
	std::ostream* report = &(std::cerr);
	inline static int nint( float in ) { return static_cast<int>(floor(in+0.5)); }

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
