#ifndef CURVE_H
#define CURVE_H

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
};


/* ----- single point data structures ----- */
//template <class T> class Curve;
//class KDeriv;
//class Parabola;
class Point {
public:
	Point( float xin=NaN, float yin=NaN, float sdensin=1. ) 
		: x(xin), y(yin), sdensity(sdensin) {}

	Point( const std::string& input ) {
		int nrd = sscanf(input.c_str(), "%f %f %f", &x, &y, &sdensity);
		if( nrd < 2 )
			throw ErrorCv::BadInput( FuncName, "format error in string "+input );
	}

	friend bool operator< ( const Point& p1, const Point& p2 ) {
		return (p1.x < p2.x);
	}

	friend std::ostream& operator<< ( std::ostream& o, const Point& pt ) {
		o << pt.x << " " << pt.y << " " << pt.sdensity;
		return o;
	}

	//friend Curve<Point>;
	//friend KDeriv;
	//friend Parabola;
	//friend Curve<Point> operator-(const Curve<Point>& c1, const Curve<Point>& c2);

	static constexpr float NaN = -12345.;
	static constexpr float twopi = M_PI * 2.;

//protected:
	float x = NaN, y = NaN, sdensity = 1.;
};

/* ----- data containers ----- */
template <class T>
class Curve {
public:
		// con/destructors
		Curve() {}
		Curve( const std::string& fname ) { Load(fname); }
		void Load( const std::string& fname ) {
			// check infile
			std::ifstream fin(fname);
			if( ! fin )
				throw ErrorCv::BadFile( FuncName, "Dispersion infile "+fname );
			// read
			for(std::string line; std::getline(fin, line); ) {
				dataV.push_back( T(line) );
			}
			std::cout<<dataV.size()<<" points loaded from file "<<fname<<std::endl;
			// sort
			Sort();
			// debug
			//for( const auto& d : dataV ) std::cerr<<d<<std::endl;
		}

		// data management and getter
		void clear() { dataV.clear();	}
		size_t size() const { return dataV.size(); }
		void reserve( size_t size ) { dataV.reserve(size);	}
		void push_back( float x, float y, float om = 1. ){ dataV.push_back( Point(x, y, om) ); }
		void Sort() { 
			std::sort(dataV.begin(), dataV.end()); 
			nsorted = dataV.size();
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
			if( itupp == dataV.end() ) {
				if( fabs((*(itupp-1)).x - x) < 0.1 ) return (*(itupp-1)).y;
				else return T::NaN;
			}
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return T::NaN;
			const auto itlow = itupp - 1;
			float vdiff = (*itupp).y - (*itlow).y, Tdiff = (*itupp).x - (*itlow).x;
			float deriv = vdiff / Tdiff;
			return (*itlow).y + ( x-(*itlow).x ) * deriv;
		}

		// find all intersections of the curve with a given horizontal line
		void FindValue(float val, std::vector<float>& perV) {
			perV.clear();
			size_t dsize = dataV.size();
			if( dsize <= 1 ) return;
			if( dsize > nsorted ) Sort();
			for( int i=0; i<dsize-1; i++ ) {
				const auto& datacur = dataV[i];
				float velh = dataV[i+1].y, vell = datacur.y;
				if( val == vell ) {	// equal to or
					perV.push_back(datacur.x);
				} else if( (val-vell) * (val-velh) < 0. ) { // in between
					float perh = dataV[i+1].x, perl = datacur.x;
					float vdiff = velh - vell, pdiff = perh - perl;
					//std::cerr<<"looking for "<<val<<" in between "<<perl<<"-"<<perh<<"sec "<<vell<<"-"<<velh<<"s/km"<<std::endl;
					perV.push_back( perl + (val-vell) * pdiff / vdiff );
				}
			}
			const auto& datacur = dataV.back();
			if( val == datacur.y ) perV.push_back(datacur.x);
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
		void Write( const std::string& fname ) {
			std::ofstream fout( fname );
			if( ! fout )
				throw ErrorCv::BadFile( FuncName, "writing to file "+fname );
			for( const auto& spP : dataV )
				fout<<spP<<"\n";
		}

		// math operators
		Curve<T>& operator-=( const Curve<T>& c2 ) { 
			(*this) = (*this) - c2; 
			return *this;
		}
		friend const Curve<T> operator-(const Curve<T>& c1, const Curve<T>& c2) {
			Curve<T> c3; c3.reserve( c1.size() );
			for( auto p1 : c1.dataV ) {
				float val2 = c2.Val(p1.x);
				if( val2 == Point::NaN ) continue;
				p1.y -= val2;
				c3.dataV.push_back(p1);
			}
			return c3;
		}

protected:
		std::vector<T> dataV;
		//std::deque<T> dataV;
private:
		typename std::vector<T>::iterator iter = dataV.begin();
		size_t nsorted = 0;
};



class Parabola {
public:
	Parabola( const Point& P1in, const Point& P2in, const Point& P3in )
		: P1(P1in), P2(P2in), P3(P3in) {}

	void Solve(){
		float x1 = P1.x, y1 = P1.y;
		float x2 = P2.x, y2 = P2.y;
		float x3 = P3.x, y3 = P3.y;
		float xs1 = x1*x1, xs2 = x2*x2, xs3 = x3*x3;
		float denom = (x1-x2) * (x1-x3) * (x2-x3);
		a = (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)) / denom;
		b = (xs3*(y1-y2) + xs2*(y3-y1) + xs1*(y2-y3)) / denom;
		c = (x2*x3*(x2-x3)*y1 + x3*x1*(x3-x1)*y2 + x1*x2*(x1-x2)*y3) / denom;
	}

	Point Vertex() {
		if(a==NaN || b==NaN || c==NaN)
			Solve();
		if( PV.x==NaN || PV.y==NaN ) {
			PV.x = - b / (2.*a);
			PV.y = c - a*PV.x*PV.y;
		}
		return PV;
	}

	float A() { 
		if(a==NaN) Solve(); 
		return a;
	}
	float B() { 
		if(b==NaN) Solve(); 
		return b;
	}
	float C() { 
		if(c==NaN) Solve(); 
		return c;
	}

private:
	Point P1, P2, P3;
	Point PV;
	float NaN = Point::NaN;
	float a=NaN, b=NaN, c=NaN;
};

#endif
