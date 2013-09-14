#include <stdio.h>
#include <iostream>
#include <limits.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <sys/time.h>
#include "VoronoiDiagramGenerator.h"
using namespace std;

uint64_t ClockGetTime()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
}

int main(int argc,char **argv) 
{	
   if(argc!=2) {
      cout<<"Usage: "<<argv[0]<<" [out Source file]"<<endl;
      exit(-1);
   }

   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::uniform_real<> uni(0,1);
   boost::variate_generator< RNGType, boost::uniform_real<> > rand(gener, uni);
   FILE *ff;
   int i, ndat = 100;
   float xdat[ndat], ydat[ndat];
   struct Point S[ndat];
/* read in from source file:
   char buff[300];
   if((ff=fopen(argv[1], "r"))==NULL) {
      cout<<"Cannot access file Plot/S.txt"<<endl;
      exit(0);
   }
   for(i=0; i<ndat && fgets(buff, 300, ff)!=NULL; i++) sscanf(buff, "%f %f", &xdat[i], &ydat[i]);
   fclose(ff);
   ndat = i;
*/
// or produce source file:
   if((ff=fopen(argv[1], "w"))==NULL) {
      cout<<"Cannot open file "<<argv[1]<<endl;
      exit(0);
   }
   for(i=0; i<ndat; i++) {
      S[i].x = rand(); xdat[i] = S[i].x;
      S[i].y = rand(); ydat[i] = S[i].y;
      fprintf(ff, "%f %f\n", S[i].x, S[i].y);
   }
   fclose(ff);
//call Voronoi digram generator:
   VoronoiDiagramGenerator vdg;
   vdg.generateVoronoi(xdat,ydat,ndat, 0,1,0,1,0.);
   vdg.resetIterator();

   struct Point S1, S2;
   cout.precision(8);
   while(vdg.getNext(S1.x,S1.y,S2.x,S2.y)) cout<<S1.x<<" "<<S1.y<<endl<<S2.x<<" "<<S2.y<<endl<<">"<<endl;

   return 0;
}



