#include "Map.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class Point3 : public Point<float> {
public:
	float z = -1.;
	Point3(const float x, const float y, const float z) 
		: Point<float>(x,y), z(z) {}
};

int main( int argc, char* argv[] ) {
   if( argc!=4 && argc!=5 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [vel map] [hdis] ([reclon] [reclat])/([receiver infile (reclon reclat)])"<<std::endl;
      exit(-1);
   }

   /* read in receiver locations */
	std::vector<Point3> reclst;
	if( argc == 4 ) {
	   std::ifstream fin(argv[3]);
		if( ! fin ) {
		   std::cerr<<"Cannot read from file "<<argv[3]<<std::endl;
			exit(0);
	   }
	   for(std::string line; std::getline(fin, line); ) {
		   float reclon, reclat, fcomp;
			if(! sscanf( line.c_str(), "%f %f %f", &reclon, &reclat, &fcomp ) ) continue;
	      if( reclon < 0. ) reclon += 360.;
		   reclst.push_back( Point3(reclon, reclat, fcomp) );
	   }
		fin.close();
	} else {
		float reclon = atof(argv[3]);
		float reclat = atof(argv[4]);
	   if( reclon < 0. ) reclon += 360.;
		reclst.push_back(Point3(reclon, reclat, -1.));
	}

   /* compute map average on each receiver location */
   Map map(argv[1]);
   float hdis = atof(argv[2]);
   if( hdis <= 0. ) exit(0);
   for(size_t i=0; i<reclst.size(); i++) {
      float res, weit;
      res = map.PointAverage(reclst[i], hdis, weit);
      std::cout<<reclst[i]<<"\t"<<res<<" "<<reclst[i].z<<"\t"<<weit<<std::endl;
   }

   return 0;
}
