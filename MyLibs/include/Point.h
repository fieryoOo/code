#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <iomanip>

template <class T>
class Point {
public:
   T lon, lat;
public:
   Point(T lonin = -12345., T latin = -12345.)
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

   friend std::ostream& operator << (std::ostream& o, Point a) {
      //o << "(" << a.lon << ", " << a.lat<< ")"; 
      o.setf(std::ios::fixed);
      o << std::left << std::setprecision(4) << a.lon << " " << a.lat;
      return o;
   }

   friend bool operator==( const Point<T>& p1, const Point<T>& p2 ) {
      return (p1.lon==p2.lon && p1.lat==p2.lat);
   }
};

#endif
