#ifndef PATHAVERAGE_H
#define PATHAVERAGE_H

#include <iostream>
#include <ostream>

int calc_dist(double lati1, double long1, double lati2, double long2, double *dist);

template < class T >
class Point {
protected:
   T lon;
   T lat;
public:
   Point(T x = 0., T y = 0.) { lon=x; lat=y; }
   T Lat() { return lat; }
   T Lon() { return lon; }
   friend std::ostream& operator << (std::ostream& o, Point a) { o << "(" << a.lon << ", " << a.lat<< ")"; return o; }
};

class Path {
protected:
   char fname[100];
   Point<double> P1, P2;
   double dist;
public:
   Path( char *inname, Point<double> Pin1, Point<double> Pin2 ) {
      sprintf(fname, "%s", inname);
      P1 = Pin1; P2 = Pin2;
      calc_dist(P1.Lat(), P1.Lon(), P2.Lat(), P2.Lon(), &dist);
   }
   double PathAverage(double lamda);
   
};


#endif
