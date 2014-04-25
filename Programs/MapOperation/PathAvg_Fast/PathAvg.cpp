#include "Map.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int main( int argc, char* argv[] ) {
   if( argc != 6 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [srclon] [srclat] [receiver infile (reclon reclat)] [vel map] [wavelength]"<<std::endl;
      exit(-1);
   }

   /* read in source location */
   Point<float> src(atof(argv[1]), atof(argv[2]));

   /* read in receiver locations */
   std::ifstream fin(argv[3]);
   if( ! fin ) {
      std::cerr<<"Cannot read from file "<<argv[3]<<std::endl;
      exit(0);
   }
   std::vector< Point<float> > reclst;
   for(std::string line; std::getline(fin, line); ) {
      float reclon, reclat;
      if(! sscanf( line.c_str(), "%f %f", &reclon, &reclat ) ) continue;
      if( reclon < 0. ) reclon += 360.;
      reclst.push_back( Point<float>(reclon, reclat) );
   }
   fin.close();

   /* compute map average along each path */
   Map map(argv[4], src);
   float lamda = atof(argv[5]);
   if( lamda < 0. ) exit(0);
   for(size_t i=0; i<reclst.size(); i++) {
      float res, perc;
      res = map.PathAverage(reclst.at(i), lamda, perc).Data();
      std::cout<<res<<" "<<perc<<std::endl;
   }
   return 0;
}
