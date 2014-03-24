#include "Map.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int main( int argc, char* argv[] ) {
   if( argc != 4 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [receiver infile (reclon reclat)] [vel map] [hdis]"<<std::endl;
      exit(-1);
   }

   /* read in receiver locations */
   std::ifstream fin(argv[1]);
   if( ! fin ) {
      std::cerr<<"Cannot read from file "<<argv[1]<<std::endl;
      exit(0);
   }
   std::vector< Point<float> > reclst;
   for(std::string line; std::getline(fin, line); ) {
      float reclon, reclat;
      sscanf( line.c_str(), "%f %f", &reclon, &reclat );
      if( reclon < 0. ) reclon += 360.;
      reclst.push_back( Point<float>(reclon, reclat) );
   }
   fin.close();

   /* compute map average on each receiver location */
   Map map(argv[2]);
   float hdis = atof(argv[3]);
   if( hdis <= 0. ) exit(0);
   for(size_t i=0; i<reclst.size(); i++) {
      float res, weit;
      res = map.PointAverage(reclst.at(i), hdis, weit);
      std::cout<<res<<" "<<weit<<std::endl;
   }
   return 0;
}
