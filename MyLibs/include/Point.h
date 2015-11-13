#ifndef POINT_H
#define POINT_H

#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

template <class T>
class Point {
public:
   T lon, lat;

public:
   //static constexpr float maxmisloc = 0.001; // allow ~0.1km mislocation
   //static constexpr float maxmislocS = maxmisloc*maxmisloc;
   static const int NaN = -12345.;

public:
   Point( T init = NaN ) 
		: lon(init), lat(init) {}
   Point(T lonin, T latin)
      : lon(lonin), lat(latin) {}
   Point( const std::string& line ) {
      LoadLine(line);
   }

	void correctLon() { if(lon<0.) lon+=360.; }

   virtual bool LoadLine( const std::string& line ) {
      //return ( sscanf(line.c_str(), "%f %f", &lon, &lat) == 2 );
		std::stringstream ss(line);
		return (ss >> lon >> lat);
   }

	bool isValid() { return (lon!=NaN && lat!=NaN); }

   inline const T& Lat() const { return lat; }
   inline       T& Lat() { return lat; }
   inline const T& Lon() const { return lon; }
   inline       T& Lon() { return lon; }

   bool IsSameLocation( const Point& p2 ) const {
      float dislon = lon-p2.lon, dislat = lat-p2.lat;
      return (dislon*dislon + dislat*dislat < 1.0e-6);
   }

	bool isWithin( const Point& BL, const Point& TR ) const {
		return lon>=BL.lon&&lon<TR.lon && lat>=BL.lat&&lat<TR.lat;
	}

	const std::string toString() { return std::to_string(lon) + " " + std::to_string(lat); }

   friend std::ostream& operator << (std::ostream& o, Point a) {
      //o << "(" << a.lon << ", " << a.lat<< ")"; 
      o.setf(std::ios::fixed);
      o << std::left << std::setprecision(4) << a.lon << " " << a.lat;
      return o;
   }

   friend bool operator==( const Point<T>& p1, const Point<T>& p2 ) {
      return (p1.lon==p2.lon && p1.lat==p2.lat);
   }

	/* math operators */

	friend bool operator<( const Point<T>& p1, const Point<T>& p2 ) {
		return p1.lon<p2.lon;
	}

	// unary negation 
	Point<T> operator-() const { return Point<T>( -lon, -lat ); }

	// addition 
	Point<T>& operator+=( const Point<T>& p2 ) {
		lon += p2.lon; lat += p2.lat;
		return *this; 
	}
	friend Point<T> operator+( const Point<T>& p1, const Point<T>& p2 ) {
		Point<T> pres = p1;
		pres += p2;
		return pres;
	}

	// subtraction
   Point<T>& operator-=( const Point<T>& p2 ) {
		(*this) += -p2;
      return *this;
   }
   friend Point<T> operator-( const Point<T>& p1, const Point<T>& p2 ) {
		Point<T> pres = p1;
      pres -= p2;
      return pres;
   }

	// multiplication (*float)
	Point<T>& operator*=( const float mul ) { 
		lon *= mul; lat *= mul;
		return *this; 
	}
	friend Point<T> operator*( const Point<T>& p1, float mul ) {
		Point<T> pres = p1;
		pres *= mul;
		return pres;
	}
	friend Point<T> operator*( float mul, const Point<T>& p1 ) {
		Point<T> pres = p1;
		pres *= mul;
		return pres;
	}

	// multiplication (*AziData)
	Point<T>& operator*=( const Point<T>& p2 ) { 
		lon *= p2.lon; lat *= p2.lat;
		return *this; 
	}
	friend Point<T> operator*( const Point<T>& p1, const Point<T>& p2 ) {
		Point<T> pres = p1;
		pres *= p2;
		return pres;
	}

	// division (*float)
	Point<T>& operator/=( const float den ) {
		float mul = 1./den;
		*this *= mul;
		return *this;
	}
	friend Point<T> operator/( const Point<T>& p1, float den ) {
		Point<T> pres = p1;
		pres /= den;
		return pres;
	}
	friend Point<T> operator/( float den, const Point<T>& p1 ) {
		Point<T> pres;
		pres.lon = den / p1.lon;
		pres.lat = den / p1.lat;
		return pres;
	}

	// division (*AziData)
	Point<T>& operator/=( const Point<T>& p2 ) { 
		lon /= p2.lon; lat /= p2.lat;
		return *this; 
	}
	friend Point<T> operator/( const Point<T>& p1, const Point<T>& p2 ) {
		Point<T> pres = p1;
		pres /= p2;
		return pres;
	}

	// square root
	Point<T> sqrt() const { return Point<T>( std::sqrt(lon), std::sqrt(lat) ); }
	friend Point<T> sqrt( const Point<T>& p1 ) { return p1.sqrt();	}

	// fabs
	Point<T> fabs() const { return Point<T>( std::fabs(lon), std::fabs(lat) ); }
	friend Point<T> fabs( const Point<T>& p1 ) { return p1.fabs();	}
};

#endif
