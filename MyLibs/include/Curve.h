#ifndef CURVE_H
#define CURVE_H

#include "Parabola.h"
#include <iostream>
#include <fstream>
#include <vector>
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
	PointC( float xin=NaN, float yin=NaN, float sdensin=1. ) 
		: x(xin), y(yin), sdensity(sdensin) {}

	PointC( const std::string& input ) {
		int nrd = sscanf(input.c_str(), "%f %f %f", &x, &y, &sdensity);
		if( nrd < 2 )
			throw ErrorCv::BadInput( FuncName, "format error in string "+input );
	}

	friend bool operator< ( const PointC& p1, const PointC& p2 ) {
		return (p1.x < p2.x);
	}

	friend std::ostream& operator<< ( std::ostream& o, const PointC& pt ) {
		o << pt.x << " " << pt.y << " " << pt.sdensity;
		return o;
	}

	//friend Curve<PointC>;
	//friend KDeriv;
	//friend Parabola;
	//friend Curve<PointC> operator-(const Curve<PointC>& c1, const Curve<PointC>& c2);

	static constexpr float NaN = -12345.;
	static constexpr float twopi = M_PI * 2.;

//protected:
	float x = NaN, y = NaN, sdensity = 1.;
};
#endif	// POINTC

/* ----- data containers ----- */
template <class T>
class Curve {
public:
		// con/destructors
		Curve( const std::string& fname = "" ) {
			// force T to be derived from PointC at compile time
			const PointC& pc = T{};
			if( !fname.empty() ) Load(fname); 
		}
		Curve( const float init ) {
			const PointC& pc = T{};
			//dataV.resize(size);
		}
		// load curve from file
		void Load( const std::string& fname ) {
			// check infile
			std::ifstream fin(fname);
			if( ! fin )
				throw ErrorCv::BadFile( FuncName, "Dispersion infile "+fname );
			// read
			for(std::string line; std::getline(fin, line); ) {
				try {	
					T t(line); 
					dataV.push_back(t);
				} catch( std::exception& e ) {}
			}
			std::cout<<dataV.size()<<" points loaded from file "<<fname<<std::endl;
			// sort
			Sort();
			// debug
			//for( const auto& d : dataV ) std::cerr<<d<<std::endl;
		}
		// convert from thickness-val data
		void FromThickData( const std::vector<T>& dataVthk, const float x0, bool addbot = true, bool addtop = true ) {
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

		// data management and getter
		void clear() { dataV.clear();	}
		size_t size() const { return dataV.size(); }
		void reserve( size_t size ) { dataV.reserve(size);	}
		//void push_back( float x, float y, float om = 1. ){ dataV.push_back( PointC(x, y, om) ); }
		void push_back( const T& t ){ dataV.push_back(t); }
		void append( const Curve<T>& c2 ) {
			auto iter1 = dataV.end()-1, iter2 = c2.dataV.begin();
			if( iter1->x == iter2->x ) {
				iter1->y = (iter1->y + iter2->y) * 0.5;
			} else {
				iter1 = dataV.end();
			}
			dataV.insert(iter1, iter2, c2.dataV.end()); 
		}
		void Sort() { 
			std::sort(dataV.begin(), dataV.end()); 
			nsorted = dataV.size();
		}
		void cut( const float xl, const float xu ) {
			Sort();
			const auto itlo = std::lower_bound( dataV.begin(), dataV.end(), T(xl) );
			dataV.erase( dataV.begin(), itlo );
			const auto itup = std::upper_bound( dataV.begin(), dataV.end(), T(xu) );
			dataV.erase( itup, dataV.end() );
			iter = dataV.begin();
		}

		typename std::vector<T>::iterator begin() { return dataV.begin(); }
		typename std::vector<T>::iterator end() { return dataV.end(); }
		typename std::vector<T>::iterator lower_bound(float x) { return std::lower_bound(dataV.begin(), dataV.end(), T(x)); }
		typename std::vector<T>::iterator upper_bound(float x) { return std::upper_bound(dataV.begin(), dataV.end(), T(x)); }
		void rewind() { iter = dataV.begin(); }
		bool get( T& Tout ) {
			if( iter >= dataV.end() ) return false;
			Tout = *iter; return true;
		}
		void next() { iter++; }

		// x-range
		float xmin() { return dataV.front().x; }
		float xmax() { return dataV.back().x; }

		// return y value at a given x
		float Val( const float x ) const {
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
			float vdiff = (*itupp).y - (*itlow).y, Tdiff = (*itupp).x - (*itlow).x;
			float deriv = vdiff / Tdiff;
			return (*itlow).y + ( x-(*itlow).x ) * deriv;
		}

		// find all intersections of the curve with a given horizontal line
		void FindValue(float val, std::vector<float>& xV) {
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

		// compute the RMS at each point within the given bin
		void BinRMS( float xhlen, Curve<T>& curveout ) {
			// check bin width
			float xrange = xmax() - xmin();
			if( xrange <= xhlen )
				throw ErrorCv::BadParam( FuncName, "x-halflength("+std::to_string(xhlen)+") <= than x-range("+std::to_string(xrange)+")." );

			// check sort
			if( dataV.size() > nsorted ) Sort();

			// sum squares of all points in range for the 1st point
			curveout = *this;
			float ssum = 0.;
			auto iterb = dataV.begin();
			auto itere = std::upper_bound( dataV.begin(), dataV.end(), xmin() + xhlen );
			auto iter = iterb;
			for( ; iter<itere; iter++ ) {
				float val = iter->y;
				ssum += val * val;
			}

			// process each point based on the ssum of last point
			auto& dVout = curveout.dataV;
			int nfreedom = iterb - itere - 1;
			for( auto iterout=dVout.begin(); iterout<dVout.end(); iterout++ ) {
				iter = dataV.begin() + ( iterout - dVout.begin() );
				auto iterbnew = std::lower_bound( iterb, iter, iter->x - xhlen );
				auto iterenew = std::upper_bound( itere, dataV.end(), iter->x + xhlen );
				// substract out points before iterbnew
				for( iter=iterb; iter<iterbnew; iter++ ) {
					float val = iter->y;
					ssum -= val * val;
				}
				// sum in points before iterenew
				for( iter=itere; iter<iterenew; iter++ ) {
					float val = iter->y;
					ssum += val * val;
				}
				iterout->y = sqrt( ssum / (iterenew - iterbnew - 1) );
				iterb = iterbnew; itere = iterenew;
			}

		}

		// compute derivative at a given x
		float Deriv( const float x ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), T(x) );
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return T::NaN;
			const auto itlow = itupp - 1;
			return ( (*itupp).y - (*itlow).y ) / ( (*itupp).x - (*itlow).x );
		}

		// IO
		void Write( const std::string& fname ) const {
			std::ofstream fout( fname );
			if( ! fout )
				throw ErrorCv::BadFile( FuncName, "writing to file "+fname );
			for( const auto& spP : dataV )
				fout<<spP<<"\n";
		}

		// math operator bases
		template<class Functor>
		Curve<T> YTransform( Functor func ) const {
			Curve<T> c2 = *this;
			//for_each( c2.dataV.begin(), c2.dataV.end(), func );
			for( auto& p : c2.dataV ) p.y = func(p.y);
			return c2;
		}

		template<class Functor> 
		friend Curve<T> Arithmetic(const Curve<T>& c1, const Curve<T>& c2, Functor func) {
			Curve<T> c3;
			//std::cerr<<"In Arithmetic: c1size="<<c1.size()<<" c2size="<<c2.size()<<std::endl;
			if( c1.dataV.empty() ) {
				c3.reserve( c2.size() );
				for( auto p2 : c2.dataV ) {
					p2.y = func( 0., p2.y );
					c3.dataV.push_back(p2);
				}
			} else {
				c3.reserve( c1.size() );
				for( auto p1 : c1.dataV ) {
					float val2 = c2.Val(p1.x);
					if( val2 == PointC::NaN ) {
						c1.Write("debug1.txt");
						c2.Write("debug2.txt");
						throw ErrorCv::OutOfRange( FuncName, "x value = "+std::to_string(p1.x) );
					}
					p1.y = func( p1.y, val2 );
					c3.dataV.push_back(p1);
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



#endif
