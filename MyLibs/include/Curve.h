#ifndef CURVE_H
#define CURVE_H

#include "Parabola.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>

#ifndef FuncName
#define FuncName __FUNCTION__
#endif

// curve exceptions
namespace ErrorCv {
	class Base : public std::runtime_error {
	public:
		Base( const std::string funcname, const std::string message )
			: std::runtime_error( "Error(" + funcname + "): " + message ) {}
	};

	class BadFile : public Base {
	public:
		BadFile( const std::string funcname, const std::string info )
			: Base( funcname, "Cannot access file ("+info+")." ) {}
	};

	class BadInput : public Base {
	public:
		BadInput( const std::string funcname, const std::string info )
			: Base( funcname, "Bad input ("+info+")." ) {}
	};

	class BadParam : public Base {
	public:
		BadParam( const std::string funcname, const std::string info )
			: Base( funcname, "Bad parameters ("+info+")." ) {}
	};

	class OutOfRange : public Base {
	public:
		OutOfRange( const std::string funcname, const std::string info )
			: Base( funcname, "Out of range ("+info+")." ) {}
	};
};


#ifndef POINTC
#define POINTC
/* ----- single point data structures ----- */
//template <class T> class Curve;
//class KDeriv;
//class Parabola;
class PointC {
public:
	PointC( float xin=NaN, float yin=NaN, float zin=1. ) 
		: x(xin), y(yin), z(zin) {}

	PointC( const std::string &line, int ix=1, int iy=2, int iz=3 ) {
		std::map<int, float*> indexM{{ix,&x}, {iy,&y}, {iz,&z}};
		int icol = 1; std::stringstream ss(line);
		for(auto iter=indexM.begin(); iter!=indexM.end(); iter++) {
			for(; icol<=iter->first; ++icol ) {
				if( ! (ss >> *(iter->second)) )
					if( iter->second == &z ) break;
					else throw std::runtime_error("Error(PointC::PointC): format error, col num out of range ("+line+")");
			}
		}
	}

	bool isValid() const { return (x!=NaN && y!=NaN); }

	// Norm
	float PointC abs() const { return sqrt(x*x + y*y); }

	// unary negation 
	PointC operator-() const { return PointC( -x, -y, -z ); }

	// addition 
	PointC& operator+=( const PointC& p2 ) {
		x += p2.x; y += p2.y; z += p2.z;
		return *this; 
	}
	friend PointC operator+( const PointC& p1, const PointC& p2 ) {
		PointC pres = p1;
		pres += p2;
		return pres;
	}

	// subtraction
   PointC& operator-=( const PointC& p2 ) {
		(*this) += -p2;
      return *this;
   }
   friend PointC operator-( const PointC& p1, const PointC& p2 ) {
		PointC pres = p1;
      pres -= p2;
      return pres;
   }

	// multiplication (*float)
	PointC& operator*=( const float mul ) { 
		x *= mul; y *= mul; z *= mul;
		return *this; 
	}
	friend PointC operator*( const PointC& p1, float mul ) {
		PointC pres = p1;
		pres *= mul;
		return pres;
	}
	friend PointC operator*( float mul, const PointC& p1 ) {
		PointC pres = p1;
		pres *= mul;
		return pres;
	}

	friend bool operator< ( const PointC& p1, const PointC& p2 ) {
		return (p1.x < p2.x);
	}

	// point - line distance
	friend float PLDist( const PointC& p, const PointC& pl1, const PointC& pl2 ) {
		const PointC pt1 = p - pl1, pt2 = pl2 - pl1;
		return fabs( pt1.x*pt2.y - pt2.x*pt1.y ) / pt2.abs();
	}

	friend std::ostream& operator<< ( std::ostream& o, const PointC& pt ) {
		o << pt.x << " " << pt.y << " " << pt.z;
		return o;
	}

	//friend Curve<PointC>;
	//friend KDeriv;
	//friend Parabola;
	//friend Curve<PointC> operator-(const Curve<PointC>& c1, const Curve<PointC>& c2);

	static constexpr float NaN = -123456.;
	static constexpr float twopi = M_PI * 2.;

//protected:
	float x = NaN, y = NaN, z = 1.;
};
#endif	// POINTC

/* ----- data containers ----- */
template <class T>
class Curve {
public:
		// con/destructors
		// load from file
		Curve( const std::string& fname = "", int ix=1, int iy=2, int iz=3 ) {
			// force T to be derived from PointC at compile time
			const PointC& pc = T{};
			if( !fname.empty() ) Load(fname, ix, iy, iz); 
		}
		// with data vector moved in
		Curve( std::vector<T>&& dataVin ) 
			: dataV(std::move(dataVin)) {
			const PointC& pc = T{};
		}
		// initialized to a single value
		Curve( const size_t size, const T Tval = T{} ) {
			//const PointC& pc = T{};
			dataV.assign(size, Tval);
		}
		// initialized to a flat line; constructor T(x, y) is assumed valid
		Curve( const size_t size, const float xmin, const float xmax, const float yval ) {
			dataV.reserve(size);
			const float step = size==1 ? 0. : (xmax-xmin) / (size-1);
			for(int i=0; i<size; i++) {
				float x = xmin + i * step;
				dataV.push_back( T(x, yval) ); 
			}
		}
		// range
		typename std::vector<T>::const_iterator begin() const { return dataV.begin(); }
		typename std::vector<T>::const_iterator end() const { return dataV.end(); }
		typename std::vector<T>::iterator begin() { return dataV.begin(); }
		typename std::vector<T>::iterator end() { return dataV.end(); }
		typename std::vector<T>::const_iterator lower_bound(float x) const { return std::lower_bound(dataV.begin(), dataV.end(), T(x)); }
		typename std::vector<T>::const_iterator upper_bound(float x) const { return std::upper_bound(dataV.begin(), dataV.end(), T(x)); }
		typename std::vector<T>::iterator lower_bound(float x) { return std::lower_bound(dataV.begin(), dataV.end(), T(x)); }
		typename std::vector<T>::iterator upper_bound(float x) { return std::upper_bound(dataV.begin(), dataV.end(), T(x)); }

		// data getter and management
		std::vector<T>& getDataV() { return dataV; }
		void clear() { dataV.clear();	}
		size_t size() const { return dataV.size(); }
		void reserve( size_t size ) { dataV.reserve(size);	}
		//void push_back( float x, float y, float om = 1. ){ dataV.push_back( PointC(x, y, om) ); }
		void push_back( const T& t ){ dataV.push_back(t); }
		void insert( typename std::vector<T>::iterator ipos, 
						 typename std::vector<T>::iterator iib, typename std::vector<T>::iterator iie ){ 
			dataV.insert(ipos, iib, iie); 
		}
		void rewind() { iter = dataV.begin(); }
		bool get( T& Tout ) {
			if( iter >= dataV.end() ) return false;
			Tout = *iter; return true;
		}
		void next() { iter++; }

		Curve<T>& cleanup();

		// x-range
		float xmin() const { return dataV.front().x; }
		float xmax() const { return dataV.back().x; }
		float ymin() const {
			float ymin = std::numeric_limits<float>::max();
			for(const auto &t : dataV) if(t.y<ymin) ymin=t.y;
			return ymin;
		}
		float ymax() const {
			float ymax = std::numeric_limits<float>::lowest();
			for(const auto &t : dataV) if(t.y>ymax) ymax=t.y;
			return ymax;
		}

		// load curve from file
		void Load( const std::string& fname, int ix=1, int iy=2, int iz=3 );

		// convert from thickness-val data
		void FromThickData( const std::vector<T>& dataVthk, const float x0, bool addbot = true, bool addtop = true );

		void append( const Curve<T>& c2 );
		void Sort();
		void cut( const float xl, const float xu );

		// return y value at a given x
		float Val( const float x, bool allowOOB = false ) const;

		// return z value at a given x
		float Z( const float x ) const;

		T GetPoint( const float x, bool allowOOB = false ) const;

		// find all intersections of the curve with a given horizontal line
		void FindValue(float val, std::vector<float>& xV);

		// curve simplification with the Ramer-Douglas-Peucker algorithm
		Curve<T> RDPSimplification( typename std::vector<T>::const_iterator iterb,
											 typename std::vector<T>::const_iterator itere, float epsilon ) const;
		Curve<T> RDPSimplification( float epsilon ) const {
			return RDPSimplification( begin(), end()-1, epsilon );
		}

		// compute the RMS at each point within the given bin
		Curve<T> BinRMS( float xhlen );

		// compute derivative/slope at a given x
		float Deriv( const float x );

		// compute slope/derivative by fitting a parabola within a given segment
		// CAUTION: dataV[].z s will be modified (for weighting!)
		float Slope( const float x, const float hlen, float& yval, float& rms, bool gaussWeight = true ) {
			return FitParabola(x, hlen, yval, rms, gaussWeight).Slope(x);
		}
		Parabola FitParabola( const float x, const float hlen, float& yval, float& rms, bool gaussWeight = true );

		// compute slope/derivative by fitting a line within a given segment
		float SlopeL( const float x, const float hlen, float& yval, float& sigma ) const;

		// compute the least-square line fit in a given range
		void FitLine(typename std::vector<T>::const_iterator iterl, 
						 typename std::vector<T>::const_iterator iteru, 
						 double &aout, double &sigma_a, double &bout, double &sigma_b ) const;

		// merge with another curve (y=current_curve, z=new_curve)
		void Merge( const Curve<T> c2 );

		// IO
		void Write( const std::string& fname ) const;


		// math operator bases
		template<class Functor>
		Curve<T> YTransform( Functor func ) const;

		template<class Functor> 
		friend Curve<T> Arithmetic(const Curve<T>& c1, const Curve<T>& c2, Functor func) {
			Curve<T> c3;
			//std::cerr<<"In Arithmetic: c1size="<<c1.size()<<" c2size="<<c2.size()<<std::endl;
			if( c1.dataV.empty() ) {
				c3.reserve( c2.size() );
				for( auto p2 : c2.dataV ) {
					p2.y = func( 0., p2.y );
					c3.push_back(p2);
				}
			} else {
				c3.reserve( c1.size() );
				for( auto p1 : c1 ) {
					float val2 = c2.Val(p1.x);
					if( val2 == PointC::NaN ) {
						//throw ErrorCv::OutOfRange( FuncName, "x value = "+std::to_string(p1.x) );
						p1.y = PointC::NaN;
					} else {
						p1.y = func( p1.y, val2 );
					}
					c3.push_back(p1);
				}
			}
			return c3;
		}


		// math operators
		Curve<T> operator-() const { return YTransform( std::negate<float>() ); }
		Curve<T> sqrt() const { return YTransform( [](const float val){ return std::sqrt(val); } ); }
		Curve<T>& operator+=( const Curve<T>& c2 ) { *this = (*this) + c2; return *this; }
		Curve<T>& operator-=( const Curve<T>& c2 ) {	*this = (*this) - c2; return *this; }
		Curve<T>& operator*=( const Curve<T>& c2 ) {	*this = (*this) * c2; return *this; }
		Curve<T>& operator/=( const Curve<T>& c2 ) {	*this = (*this) / c2; return *this; }
		Curve<T>& operator*=( const float mul ) { 
			for( auto& p : dataV ) p.y *= mul;
			return *this;
		}
		Curve<T>& operator/=( const float div ) {
			return ( (*this) *= (1./div) );
		}

		friend Curve<T> sqrt( const Curve<T>& c ) { return c.sqrt(); }
		friend Curve<T> operator+(const Curve<T>& c1, const Curve<T>& c2) {
			return Arithmetic( c1, c2, std::plus<float>() );
		}
		friend Curve<T> operator-(const Curve<T>& c1, const Curve<T>& c2) {
			return Arithmetic( c1, c2, std::minus<float>() );
		}
		friend Curve<T> operator*(const Curve<T>& c1, const Curve<T>& c2) {
			return Arithmetic( c1, c2, std::multiplies<float>() );
		}
		friend Curve<T> operator*(const Curve<T>& c1, const float mul) {
			return c1.YTransform( [&](const float val){ return val*mul; } );
		}
		friend Curve<T> operator/(const Curve<T>& c1, const Curve<T>& c2) {
			return Arithmetic( c1, c2, std::divides<float>() );
		}
		friend Curve<T> operator/(const Curve<T>& c1, const float div) {
			return c1 * (1./div);
		}

protected:
		std::vector<T> dataV;
		//std::deque<T> dataV;
private:
		typename std::vector<T>::iterator iter = dataV.begin();
		size_t nsorted = 0;
};


/* -------------------- Implementations -------------------- */

// load curve from file
template <class T>
void Curve<T>::Load( const std::string& fname, int ix, int iy, int iz ) {
	// check infile
	std::ifstream fin(fname);
	if( ! fin )
		throw ErrorCv::BadFile( FuncName, "Dispersion infile "+fname );
	// read
	for(std::string line; std::getline(fin, line); ) {
		try {	
			T t(line, ix, iy, iz); 
			dataV.push_back(t);
		} catch( std::exception& e ) {}
	}
	//std::cout<<dataV.size()<<" points loaded from file "<<fname<<std::endl;
	// sort
	Sort();
	// debug
	//for( const auto& d : dataV ) std::cerr<<d<<std::endl;
}

// convert from thickness-val data
template <class T>
void Curve<T>::FromThickData( const std::vector<T>& dataVthk, const float x0, bool addbot, bool addtop ) {
	const size_t thksize = dataVthk.size();
	//std::cerr<<"bb thksize = "<<thksize<<"  datasize = "<<dataV.size()<<std::endl;
	const size_t newsize = thksize + addtop + addbot;
	dataV.resize(newsize);
	nsorted = 0; iter = dataV.begin();
	// current depth
	float dep = x0;
	// first (bottom) point
	if( addbot ) {
		dataV[0] = dataVthk[0]; dataV[0].x = dep;
	}
	// all thk points
	for( int i=0; i<thksize; i++ ) {
		auto& p = dataV[i+addbot];
		p = dataVthk[i];
		float half_thk = p.x*0.5;
		dep += half_thk;
		p.x = dep;
		dep += half_thk;
	}
	// last (top) point
	if( addtop ) {
		dataV[newsize-1] = dataVthk[thksize-1]; dataV[newsize-1].x = dep;
	}
}

template <class T>
void Curve<T>::append( const Curve<T>& c2 ) {
	auto iter1 = dataV.end()-1, iter2 = c2.dataV.begin();
	if( iter1->x == iter2->x ) {
		iter1->y = (iter1->y + iter2->y) * 0.5;
	} else {
		iter1 = dataV.end();
	}
	dataV.insert(iter1, iter2, c2.dataV.end()); 
}
template <class T>
void Curve<T>::Sort() { 
	std::sort(dataV.begin(), dataV.end()); 
	nsorted = dataV.size();
}
template <class T>
void Curve<T>::cut( const float xl, const float xu ) {
	Sort();
	const auto itlo = std::lower_bound( dataV.begin(), dataV.end(), T(xl) );
	dataV.erase( dataV.begin(), itlo );
	const auto itup = std::upper_bound( dataV.begin(), dataV.end(), T(xu) );
	dataV.erase( itup, dataV.end() );
	iter = dataV.begin();
}

// return y value at a given x
template <class T>
float Curve<T>::Val( const float x, bool allowOOB ) const {
	const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), T(x) );
	// greater than end
	const float ferr = 1.0e-3;
	if( itupp >= dataV.end() ) {
		if( allowOOB || fabs((*(itupp-1)).x - x)<ferr ) return (*(itupp-1)).y;
		else return T::NaN;
	}
	// smaller than beg
	if( itupp<dataV.begin()+1 ) {
		if( allowOOB || fabs((*(itupp)).x - x)<ferr ) return (*(itupp)).y;
		return T::NaN;
	}
	const auto itlow = itupp - 1;
	if( itlow->x > x+ferr )	return T::NaN;
	// between two points: interpolate
	float vdiff = (*itupp).y - (*itlow).y, Tdiff = (*itupp).x - (*itlow).x;
	float deriv = vdiff / Tdiff;
	return (*itlow).y + ( x-(*itlow).x ) * deriv;
}

// return z value at a given x
template <class T>
float Curve<T>::Z( const float x ) const {
	const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), T(x) );
	// at the end
	const float ferr = 1.0e-3;
	if( itupp == dataV.end() ) {
		if( fabs((*(itupp-1)).x - x) < ferr ) return (*(itupp-1)).y;
		else return T::NaN;
	}
	// out of range
	if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
		return T::NaN;
	const auto itlow = itupp - 1;
	if( itlow->x > x+ferr )	return T::NaN;
	// between two points: interpolate
	float vdiff = (*itupp).z - (*itlow).z, Tdiff = (*itupp).x - (*itlow).x;
	float deriv = vdiff / Tdiff;
	return (*itlow).z + ( x-(*itlow).x ) * deriv;
}

template <class T>
T Curve<T>::GetPoint( const float x, bool allowOOB ) const {
	const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), T(x) );
	// x out of range:
	const float ferr = 1.0e-3;
	if( itupp == dataV.end() ) {
		if( allowOOB || fabs((*(itupp-1)).x - x)<ferr ) return (*(itupp-1));
		else throw ErrorCv::BadParam( FuncName, "requested x out of range!" );
	}
	if( itupp<dataV.begin()+1 ) {
		if( allowOOB || fabs((*(itupp)).x - x)<ferr ) return (*(itupp));
		else throw ErrorCv::BadParam( FuncName, "requested x out of range!" );
	}
	// between two points: interpolate
	const auto itlow = itupp - 1;
	return ( *itlow + (*itupp-*itlow) * (x-itlow->x) * (1./(itupp->x-itlow->x)) );
}

// find all intersections of the curve with a given horizontal line
template <class T>
void Curve<T>::FindValue(float val, std::vector<float>& xV) {
	xV.clear();
	size_t dsize = dataV.size();
	if( dsize <= 1 ) return;
	if( dsize > nsorted ) Sort();
	for( int i=0; i<dsize-1; i++ ) {
		const auto& datacur = dataV[i];
		float velh = dataV[i+1].y, vell = datacur.y;
		if( val == vell ) {	// equal to or
			xV.push_back(datacur.x);
		} else if( (val-vell) * (val-velh) < 0. ) { // in between
			float perh = dataV[i+1].x, perl = datacur.x;
			float vdiff = velh - vell, pdiff = perh - perl;
			//std::cerr<<"looking for "<<val<<" in between "<<perl<<"-"<<perh<<"sec "<<vell<<"-"<<velh<<"s/km"<<std::endl;
			xV.push_back( perl + (val-vell) * pdiff / vdiff );
		}
	}
	const auto& datacur = dataV.back();
	if( val == datacur.y ) xV.push_back(datacur.x);
}

template <class T>
Curve<T>& Curve<T>::cleanup() {
	Curve<T> c2; c2.reserve( size() );
	for( const auto& p : dataV )
		if( p.y != T::NaN ) c2.push_back(p);
	*this = std::move(c2);
	return *this;
}

template <class T>
Curve<T> Curve<T>::RDPSimplification( typename std::vector<T>::const_iterator iterb,
												  typename std::vector<T>::const_iterator itere, float epsilon ) const {
	// find the point with max distance
	auto itermax = iterb;
	float dismax = 0.;
	for(auto iter=iterb+1; iter<itere; iter++) {
		float dis = PLDist(*iter, *iterb, *itere);
		if( dis > dismax ) {
			itermax = iter; dismax = dis;
		}
	}
	//std::cerr<<*itermax<<"  -  "<<*iterb<<" "<<*itere<<std::endl;
	//std::cerr<<dismax<<std::endl;

	Curve<T> cvres;
	// recurse when dismax > epsilon
	if( dismax > epsilon ) {
		// recurse
		Curve<T> cv1 = RDPSimplification( iterb, itermax, epsilon );
		Curve<T> cv2 = RDPSimplification( itermax, itere, epsilon );
		// and merge
		cv1.reserve( cv1.size() + cv2.size() - 1 );
		std::move( cv2.begin()+1, cv2.end(), std::back_inserter(cv1.dataV) );
		cvres.dataV = std::move(cv1.dataV);
	} else {
		cvres.push_back(*iterb);
		cvres.push_back(*itere);
	}

	return cvres;
}

// compute the RMS at each point within the given bin
template <class T>
Curve<T> Curve<T>::BinRMS( float xhlen ) {
	// check bin width
	float xrange = xmax() - xmin();
	if( xrange <= xhlen )
		throw ErrorCv::BadParam( FuncName, "x-halflength("+std::to_string(xhlen)+") <= than x-range("+std::to_string(xrange)+")." );

	// check sort
	if( dataV.size() > nsorted ) Sort();

	// sum squares of all points in range for the 1st point
	Curve<T> curveout(*this);
	// be careful with rounding error!!!
	long double ssum = 0.;
	auto iterb = dataV.begin();
	auto itere = std::upper_bound( dataV.begin(), dataV.end(), xmin() + xhlen );
	auto iter = iterb;
	for( ; iter<itere; iter++ ) {
		long double val = iter->y;
		ssum += val * val;
	}

	// process each point based on the ssum of last point
	auto& dVout = curveout.dataV;
	for( auto iterout=dVout.begin(); iterout<dVout.end(); iterout++ ) {
		iter = dataV.begin() + ( iterout - dVout.begin() );
		auto iterbnew = std::lower_bound( iterb, iter, iter->x - xhlen );
		auto iterenew = std::upper_bound( itere, dataV.end(), iter->x + xhlen );
		// substract out points before iterbnew
		for( iter=iterb; iter<iterbnew; iter++ ) {
			long double val = iter->y;
			ssum -= val * val;
		}
		// sum in points before iterenew
		for( iter=itere; iter<iterenew; iter++ ) {
			long double val = iter->y;
			ssum += val * val;
		}
		int nfree = iterenew - iterbnew - 1;
		iterout->y = nfree>0 ? std::sqrt(ssum/nfree) : T::NaN;
		iterb = iterbnew; itere = iterenew;
	}

	return curveout;
}

// compute derivative/slope at a given x
template <class T>
float Curve<T>::Deriv( const float x ) {
	const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), T(x) );
	if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
		return T::NaN;
	const auto itlow = itupp - 1;
	return ( (*itupp).y - (*itlow).y ) / ( (*itupp).x - (*itlow).x );
}

// compute slope/derivative by fitting a parabola within a given segment
// CAUTION: dataV[].z s will be modified (for weighting!)
template <class T>
Parabola Curve<T>::FitParabola( const float x, const float hlen, float& yval, float& rms, bool gaussWeight ) {
	// locate boundaries
	float hlen_fit = gaussWeight ? 3.*hlen : hlen;
	const auto itlow = std::lower_bound( dataV.begin(), dataV.end(), T(x-hlen_fit) );
	const auto itupp = std::upper_bound( itlow, dataV.end(), T(x+hlen_fit) );
	//const float minhlen = hlen*0.25;
	const float minhlen = hlen*0.01;
	if( itlow->x>x-minhlen || (itupp-1)->x<x+minhlen || itupp-itlow<3 ) {
		rms = T::NaN; return Parabola();
	}
	// assign weight to dataV[].z
	if( gaussWeight ) {
		const float alpha = -0.5 / (hlen*hlen);
		for( auto iter=itlow; iter!=itupp; iter++ ) {
			float xdiff = iter->x - x;
			iter->z = exp( alpha * xdiff * xdiff );
		}
	} else {
		// assumes iter->zs are all 1 already
		//for( auto iter=itlow; iter!=itupp; iter++ ) iter->z = 1.0;
	}
	// fit the parabola
	//std::cerr<<*itlow<<" "<<*(itupp-1)<<" "<<itupp-itlow<<std::endl;
	Parabola parab; parab.Fit(itlow, itupp, rms);
	yval = parab(x);
	return parab;
}

// compute slope/derivative by fitting a line within a given segment
template <class T>
float Curve<T>::SlopeL( const float x, const float hlen, float& yval, float& sigma ) const {
	const auto itlow = std::lower_bound( dataV.begin(), dataV.end(), T(x-hlen) );
	const auto itupp = std::upper_bound( itlow, dataV.end(), T(x+hlen) );
	const float minhlen = hlen*0.01;
	//const float minhlen = hlen*0.25;
	if( itlow==dataV.end() || itlow->x>x-minhlen || (itupp-1)->x<x+minhlen ) {
		sigma = T::NaN; return T::NaN;
	}
	double a, b, sigma_a, sigma_b;
	FitLine(itlow, itupp, a, sigma_a, b, sigma_b);
	yval = a * x + b;	sigma = sigma_a; 
	return a;
}

// compute the least-square line fit in a given range
template <class T>
void Curve<T>::FitLine(typename std::vector<T>::const_iterator iterl, 
							  typename std::vector<T>::const_iterator iteru, 
							  double &aout, double &sigma_a, double &bout, double &sigma_b ) const {
	const int ndat = iteru - iterl;
	if( ndat < 2 ) {
		aout = bout = sigma_a = sigma_b = T::NaN;
		return;
	}

	int i;
	double W=0, WX=0, WY=0, WX2=0, WY2=0, WXY=0;
	for(auto iter=iterl;iter<iteru;iter++) {
		double x = iter->x, y = iter->y;
		double w = 1.;	// 1./( sigma[i] * sigma[i] );
		W += w;
		WX += w * x; WY += w * y;
		WX2 += w * x * x; WY2 += w * y * y;
		WXY += w * x * y; 
	}

	// determine a b assuming x to be independent
	double w = 1./(W*WX2-WX*WX);
	aout = (W*WXY-WX*WY) * w;
	bout = (-WX*WXY+WX2*WY) * w;
	/* // determine a b assuming 7 to be independent
		double w = 1./(W*WXY-WX*WY);
		aout=(W*WY2-WY*WY) * w;
		bout=(WY*WXY-WY2*WX) * w; */

	if( ndat == 2 ) {
		//sigma_a = sigma_b = std::numeric_limits<double>::max();
		sigma_a = sigma_b = T::NaN;
		return;
	}

	// compute uncertainty in a and b
	// define k according to ndat from t distribution
	// assuming a 95% conf is equivalent to 2 sigma conf
	double k;
	if( ndat == 3 ) k = 12.706;
	else if ( ndat == 4 ) k = 4.303;
	else if ( ndat == 5 ) k = 3.182;
	else if ( ndat == 6 ) k = 2.776;
	else k = 1.960+13.8/pow((double)ndat, 1.6); // this could be made better
	k *= 0.5; // now this is 1 sigma
	// compute uncertainties
	double S2 = 0., x_avg = WX/W, varX = 0.; 
	for(auto iter=iterl;iter<iteru;iter++) {
		double dx = iter->x - x_avg;
		double dy = iter->y - aout*iter->x - bout;
		S2 += dy * dy; varX += dx * dx;
	}
	S2 /= ndat-2.;
	// determine rms assuming x to be independent
	//sigma_a = std::sqrt(S2 / varX);
	sigma_a = std::sqrt(S2 * W * w * W) * k;		// note!!! output is standard-dev
	sigma_b = std::sqrt(S2 * WX2 * w * W) * k;	// due to the extra sqrt(*W) !
	/* // determine rms assuming y to be independent
		S2 /= a*a;
	 *sigmaaout = k * sqrt(S2 * W * w);
	 *sigmabout = k * sqrt(S2 * WY2 * w); */
}

// merge with another curve (y=current_curve, z=new_curve)
template <class T>
void Curve<T>::Merge( const Curve<T> c2 ) {
	// lambda function to insert points in the range [iterl, iteru] into dataV
	auto InsertPoints = [&](typename std::vector<T>::const_iterator iterl, 
			typename std::vector<T>::const_iterator iteru) {
		for( auto iter=iterl; iter<iteru; iter++ ) {
			T p; p.x = iter->x; p.z = iter->y;
			p.y = Val(p.x, true);
			dataV.push_back( p );
		}
	};
	// insert points with x smaller than this.xmin
	auto iter2 = std::lower_bound( c2.begin(), c2.end(), T(xmin()) );
	InsertPoints( c2.begin(), iter2 );
	// for exist points (in between this->begin() and this->end()), write c2.y into p.z
	for( auto& p : dataV ) p.z = c2.Val(p.x, true);
	// insert points with x greater than this.xmax
	iter2 = std::upper_bound( iter2, c2.end(), T(xmax()) );
	InsertPoints( iter2, c2.end() );
	// sort
	Sort();
}

// IO
template <class T>
void Curve<T>::Write( const std::string& fname ) const {
	std::ofstream fout( fname );
	if( ! fout )
		throw ErrorCv::BadFile( FuncName, "writing to file "+fname );
	for( const auto& spP : dataV )
		fout<<spP<<"\n";
}


// math operator bases
template<class T> template<class Functor>
Curve<T> Curve<T>::YTransform( Functor func ) const {
	Curve<T> c2 = *this;
	//for_each( c2.dataV.begin(), c2.dataV.end(), func );
	for( auto& p : c2.dataV ) p.y = func(p.y);
	return c2;
}


#endif
