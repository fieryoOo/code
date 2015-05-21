#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <sstream>
#include <iomanip>

template <class T>
class Point {
public:
   T lon = NaN, lat = NaN;

public:
   static constexpr float maxmisloc = 0.001; // allow ~0.1km mislocation
   static constexpr float maxmislocS = maxmisloc*maxmisloc;
   static constexpr float NaN = -12345.;

public:
   Point(T lonin = NaN, T latin = NaN)
      : lon(lonin), lat(latin) {}
   Point( const std::string& line ) {
      LoadLine(line);
   }

   virtual bool LoadLine( const std::string& line ) {
      return ( sscanf(line.c_str(), "%f %f", &lon, &lat) == 2 );
   }

   inline const T& Lat() const { return lat; }
   inline       T& Lat() { return lat; }
   inline const T& Lon() const { return lon; }
   inline       T& Lon() { return lon; }

   bool IsSameLocation( const Point& p2 ) const {
      float dislon = lon-p2.lon, dislat = lat-p2.lat;
      return (dislon*dislon + dislat*dislat < maxmislocS);
   }

   friend std::ostream& operator << (std::ostream& o, Point a) {
      //o << "(" << a.lon << ", " << a.lat<< ")"; 
      o.setf(std::ios::fixed);
      o << std::left << std::setprecision(4) << a.lon << " " << a.lat;
      return o;
   }

   friend bool operator==( const Point<T>& p1, const Point<T>& p2 ) {
      return (p1.lon==p2.lon && p1.lat==p2.lat);
   }

	friend bool operator<( const Point<T>& p1, const Point<T>& p2 ) {
		return p1.lon<p2.lon;
	}
};

#endif
