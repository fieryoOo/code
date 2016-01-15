#ifndef PARABOLA_H
#define PARABOLA_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdexcept>

#ifndef FuncName
#define FuncName __FUNCTION__
#endif

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

	PointC( const std::string& input ) {
		int nrd = sscanf(input.c_str(), "%f %f %f", &x, &y, &z);
		if( nrd < 2 )
			throw std::runtime_error( std::string(FuncName) + ": format error in string "+input );
	}

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

	friend std::ostream& operator<< ( std::ostream& o, const PointC& pt ) {
		o << pt.x << " " << pt.y << " " << pt.z;
		return o;
	}

	//friend Curve<PointC>;
	//friend KDeriv;
	//friend Parabola;
	//friend Curve<PointC> operator-(const Curve<PointC>& c1, const Curve<PointC>& c2);

	static constexpr float NaN = -12345.;
	static constexpr float twopi = M_PI * 2.;

//protected:
	float x = NaN, y = NaN, z = 1.;
};
#endif	// POINTC


// a parabola of the standard form y = ax**2 + bx + c
class Parabola {
public:
	Parabola( const PointC& P1in, const PointC& P2in, const PointC& P3in )
		: P1(P1in), P2(P2in), P3(P3in) {}

	void Solve() const{
		long double x1 = P1.x, y1 = P1.y;
		long double x2 = P2.x, y2 = P2.y;
		long double x3 = P3.x, y3 = P3.y;
		long double xs1 = x1*x1, xs2 = x2*x2, xs3 = x3*x3;
		long double denom = (x1-x2) * (x1-x3) * (x2-x3);
		a = (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)) / denom;
		b = (xs3*(y1-y2) + xs2*(y3-y1) + xs1*(y2-y3)) / denom;
		c = (x2*x3*(x2-x3)*y1 + x3*x1*(x3-x1)*y2 + x1*x2*(x1-x2)*y3) / denom;
		//std::cerr<<"results: "<<a<<" "<<b<<" "<<c<<std::endl;
	}

	PointC Vertex() const {
		if(a==NaN || b==NaN || c==NaN) Solve();
		if( PV.x==NaN || PV.y==NaN ) {
			PV.x = - b / (2.*a);
			PV.y = c - a*PV.x*PV.x;
		}
		return PV;
	}

	long double A() const { 
		if(a==NaN) Solve(); 
		return a;
	}
	long double B() const { 
		if(b==NaN) Solve(); 
		return b;
	}
	long double C() const { 
		if(c==NaN) Solve(); 
		return c;
	}

	long double operator()( const long double x ) const {
		if(a==NaN || b==NaN || c==NaN) Solve();
		return ( a*x*x + b*x + c );
	}

protected:
	static constexpr long double NaN = PointC::NaN;

private:
	const PointC P1, P2, P3;
	mutable PointC PV;
	mutable long double a=NaN, b=NaN, c=NaN;
};

#endif
