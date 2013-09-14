#include <stdio.h>
#include <iostream>
#include <limits.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <sys/time.h>
#include "VoronoiFortunes.h"
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
      cout<<"Usage: "<<argv[0]<<" [out_file for time space behaviors]"<<endl;
      exit(-1);
   }

   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::uniform_real<> uni(0,1);
   boost::variate_generator< RNGType, boost::uniform_real<> > rand(gener, uni);
   FILE *ff;
   int i, j, k, ndat;
/* read in from source file:
   char buff[300];
   if((ff=fopen(argv[1], "r"))==NULL) {
      cout<<"Cannot access file Plot/S.txt"<<endl;
      exit(0);
   }
   for(i=0; i<ndat && fgets(buff, 300, ff)!=NULL; i++) sscanf(buff, "%lf %lf", &xdat[i], &ydat[i]);
   fclose(ff);
   ndat = i;
//
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
*/
//call Voronoi digram generator:
   double *xdat=NULL, *ydat=NULL;
   VoronoiFortunes *vdg = NULL;
   int niter = 3;
   float Onum, OnumVL, OnumBL, Msize, MsizeVL, MsizeBL;
//   vdg.VoronoiGenerator(xdat,ydat,ndat, 0,1,0,1);
   FILE *fout = fopen(argv[1], "a");
   for(i=1; i<15; i++) {
      ndat = 1<<i ;
      xdat = new double[ndat];
      ydat = new double[ndat];
      
      Onum = 0, OnumVL = 0, OnumBL = 0, Msize = 0, MsizeVL = 0, MsizeBL = 0;
      for(k=0; k<niter; k++) {
cout<<ndat<<" points. the "<<k<<"th iteration"<<endl;
         for(j=0; j<ndat; j++) {
            xdat[j] = rand();
            ydat[j] = rand()*100.;
         }
	 vdg = new VoronoiFortunes;
         vdg->VoronoiGenerator(xdat,ydat,ndat, 0,1.,0,100);
	 Onum += vdg->Onum; OnumVL += vdg->OnumVL; OnumBL += vdg->OnumBL;
	 Msize += vdg->Msize; MsizeVL += vdg->MsizeVL; MsizeBL += vdg->MsizeBL;
         delete vdg; vdg = NULL;
      }
      delete[] xdat; xdat = NULL;
      delete[] ydat; ydat = NULL;
      Onum /= niter; OnumVL /= niter; OnumBL /= niter;
      Msize /= niter; MsizeVL /= niter; MsizeBL /= niter;
      fprintf(fout, "%d %.1f %.1f %.1f %.1f %.1f %.1f\n", ndat, Onum, OnumVL, OnumBL, Msize, MsizeVL, MsizeBL);
   }
   fclose(fout);

//   struct Point S1, S2;
//   vdg.Rewind();
//   cout.precision(8);
//   while(vdg.GetBisector(&S1.x,&S1.y,&S2.x,&S2.y)) cout<<S1.x<<" "<<S1.y<<endl<<S2.x<<" "<<S2.y<<endl<<">"<<endl;

   return 1;
}



