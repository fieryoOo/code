#ifndef PATHAVERAGE_H
#define PATHAVERAGE_H

#include "DisAzi.h"
#include <iostream>
#include <ostream>

template < class T >
class Point {
protected:
   T lon;
   T lat;
public:
   Point(T x = 0., T y = 0.) { lon=x; lat=y; }
   T Lat() { return lat; }
   T Lon() { return lon; }
   void SetLat(T latin) { lat = latin; }
   void SetLon(T lonin) { lon = lonin; }
   void move(T dlon, T dlat) { lon += dlon; lat += dlat; }
   friend std::ostream& operator << (std::ostream& o, Point a) { o << a.lon << " " << a.lat; return o; }
};

class Path {
protected:
   char fname[100];
   Point<float> P1, P2;
   double dist;
public:
   Path( char *inname = NULL, Point<float> Pin1 = Point<float>(0., 0.), Point<float> Pin2 = Point<float>(0., 0.) ) {
      if(inname) sprintf(fname, "%s", inname);
      P1 = Pin1; P2 = Pin2;
      calc_dist((double)P1.Lat(), (double)P1.Lon(), (double)P2.Lat(), (double)P2.Lon(), &dist);
   }
   double Dist(){ return dist; }
   double PathAverage(double lamda);
   
};


#endif
