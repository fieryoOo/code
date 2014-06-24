/* this code takes (1) a list of SAC files (just to be convenient. The 2 station locations stored in the SAC header will be used to predict the arrival time) 
 * (2) a group velocity map of the region of interests, and (3) a source location
 * and predict the arrival time of the precursoring signal for each path */
#include "SACREC.h"
#include "PathAverage.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <new>
#include <vector>

int main( int argc, char* argv[] ) 
{
   /* check for and read in input parameters */
   if( argc != 4 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [SAC list (fsac source_lon source_lat)] [Input model] [uncertainty in group speed map (0. - 1.)]"<<std::endl;
      exit(-1);
   }
   char mapname[100]; sprintf(mapname, "%s", argv[2]);
   float unc = atof(argv[3]);

   /* read in sac list */
   std::vector<std::string> filebuff, namebuff;
   std::ifstream fin(argv[1]);
   if( !fin ) {
      std::cerr<<" Cannot access sac list file "<<argv[1]<<std::endl;
      exit(0);
   }
   char buff[300];
   while( fin.getline(buff, 300) ) filebuff.push_back(buff);
   fin.close();

   /* process each sac file */
   SACREC *sigcur = NULL;
   Path *pathcur = NULL;
   std::vector<std::string>::iterator ifile;
   std::vector< std::vector<float> > arrtime;
   arrtime.reserve(filebuff.size());
   namebuff.reserve(filebuff.size());
   std::vector<float> timetmp(5, 0.);
// source point
   for(ifile=filebuff.begin(); ifile<filebuff.end(); ifile++) {
      /* get fsac name and source location */
      float slon, slat;
      if( (sscanf( (*ifile).c_str(), "%s %f %f", buff, &slon, &slat)) != 3 ) continue;
      namebuff.push_back(buff);
      Point<float> Psource(slon, slat);
      /* initialize sac record and load sac header&signal to get station locations */
      sigcur = new SACREC(buff, 0., 0.);
      if( ! sigcur->load() ) {
         delete sigcur; sigcur = NULL;
         std::cerr<<"Warning(main): sac file not exist OR bad data in sac file "<<*ifile<<std::endl;
         continue;
      }
      /* compute average velocity from source to sta1 */
      pathcur = new Path(mapname, Psource, sigcur->P1());
      double t1 = pathcur->Dist() / pathcur->PathAverage(20.);
      delete pathcur; pathcur = NULL;
      /* compute average velocity from source to sta2 */
      pathcur = new Path(mapname, Psource, sigcur->P2());
      double t2 = pathcur->Dist() / pathcur->PathAverage(20.);
      delete pathcur; pathcur = NULL;
      /* compute a time range of the precursor arrival and store */
      timetmp[0] = t2*(1.+unc) - t1*(1.-unc);
      timetmp[1] = t2*(1.-unc) - t1*(1.+unc);
      /* compute average velocity from sta1 to sta2 */
      pathcur = new Path(mapname, sigcur->P1(), sigcur->P2());
      /* save distance */
      timetmp[4] = pathcur->Dist();
      double t0 = timetmp[4] / pathcur->PathAverage(20.);
      delete pathcur; pathcur = NULL;
      /* compute a time range of the surface wave arrival and store */
      timetmp[2] = t0*(1.+unc);
      timetmp[3] = t0*(1.-unc);

      arrtime.push_back(timetmp);

      delete sigcur; sigcur = NULL;
   }

   /* output results */
   char outname[150];
   sprintf(outname, "%s_Tprecursor", argv[1]);
   std::ofstream fout(outname);
   for(int i=0; i<namebuff.size(); i++) fout<<namebuff.at(i)<<" "<<arrtime.at(i).at(0)<<" "<<arrtime.at(i).at(1)<<" "<<arrtime.at(i).at(2)<<" "<<arrtime.at(i).at(3)<<" "<<arrtime.at(i).at(4)<<std::endl;
   fout.close();

   return 0;
}
