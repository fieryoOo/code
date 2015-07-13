#ifndef PARABOLA_H
#define PARABOLA_H

#include <iostream>
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
	PointC( float xin=NaN, float yin=NaN, float sdensin=1. ) 
		: x(xin), y(yin), sdensity(sdensin) {}

	PointC( const std::string& input ) {
		int nrd = sscanf(input.c_str(), "%f %f %f", &x, &y, &sdensity);
		if( nrd < 2 )
			throw std::runtime_error( std::string("Error(") + FuncName + "): Bad input (format error in string " + input + ")." );
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


// a parabola of the standard form y = ax**2 + bx + c
class Parabola {
public:
	Parabola( const PointC& P1in, const PointC& P2in, const PointC& P3in )
		: P1(P1in), P2(P2in), P3(P3in) {}

	void Solve() const{
		float x1 = P1.x, y1 = P1.y;
		float x2 = P2.x, y2 = P2.y;
		float x3 = P3.x, y3 = P3.y;
		float xs1 = x1*x1, xs2 = x2*x2, xs3 = x3*x3;
		float denom = (x1-x2) * (x1-x3) * (x2-x3);
		a = (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)) / denom;
		b = (xs3*(y1-y2) + xs2*(y3-y1) + xs1*(y2-y3)) / denom;
		c = (x2*x3*(x2-x3)*y1 + x3*x1*(x3-x1)*y2 + x1*x2*(x1-x2)*y3) / denom;
	}

	PointC Vertex() const {
		if(a==NaN || b==NaN || c==NaN) Solve();
		if( PV.x==NaN || PV.y==NaN ) {
			PV.x = - b / (2.*a);
			PV.y = c - a*PV.x*PV.y;
		}
		return PV;
	}

	float A() const { 
		if(a==NaN) Solve(); 
		return a;
	}
	float B() const { 
		if(b==NaN) Solve(); 
		return b;
	}
	float C() const { 
		if(c==NaN) Solve(); 
		return c;
	}

	float operator()( const float x ) const {
		if(a==NaN || b==NaN || c==NaN) Solve();
		return ( a*x*x + b*x + c );
	}

protected:
	static constexpr float NaN = PointC::NaN;

private:
	PointC P1, P2, P3;
	mutable PointC PV;
	mutable float a=NaN, b=NaN, c=NaN;
};

#endif
