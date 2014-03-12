#ifndef MAP_H
#define MAP_H

#include "Array2D.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

template < class T >
class Point {
   T lon, lat;
public:
   Point(T lonin = -12345., T latin = -12345.) 
      : lon(lonin), lat(latin) {}
   inline const T Lat() const { return lat; }
   inline	T Lat() { return lat; }
   inline const T Lon() const { return lon; }
   inline	T Lon() { return lon; }
   friend std::ostream& operator << (std::ostream& o, Point a) { 
      o << "(" << a.lon << ", " << a.lat<< ")"; 
      return o; 
   }
};

#define PIO180 0.017453292519943295
template < class T >
class DataPoint {
   T dis, disa, azi, data;
public:
   DataPoint(const T disin = -12345., const T aziin = -12345., const T datain = -12345.)
      : dis(disin), azi(aziin), data(datain) {
      disa = dis * azi * PIO180;
   }

   inline const T Dis() const { return dis; }
   inline	T Dis() { return dis; }
   inline const T Disa() const { return disa; }
   inline	T Disa() { return disa; }
   inline const T Azi() const { return azi; }
   inline	T Azi() { return azi; }
   inline const T Data() const { return data; }
   inline	T Data() { return data; }

   friend std::ostream& operator << (std::ostream& o, DataPoint a) { 
      o << "(" << a.dis << ", " << a.disa<< "): "<<a.data; 
      return o; 
   }
};


class Map {
   struct Mimpl;
   std::unique_ptr<Mimpl> pimplM;
/*
   std::string fname;
   Point<float> src;
   float dismax, disamax;
   std::vector< DataPoint<float> > dataV;
   Array2D< DataPoint<float> > dataM;
*/
public:
   Map( const char *inname, const Point<float>& srcin );
   ~Map();

   double PathAverage(Point<float> Prec, double lamda);
   
};


#endif
