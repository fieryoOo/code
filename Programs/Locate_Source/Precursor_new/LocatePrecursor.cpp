/* This code takes (1) a list of sac envelope files along with their center period (and optionally the measured group vel) and (2) an input map
 * of group vel model for the region of interests and computes a map of probability of the precursory signal location. */
#include "PathAverage.h"
#include "SACREC.h"
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <new>
#include <vector>

#define MAIN
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_num_threads() { return 1;}
#endif

/* input parameters:
 * 1) Input list: col1 = sac envelope list, col2 = center period, col3(optional) = measured group vel of this path
 * 2) Group vel map: col1 = lon, col2 = lat, col3 = grvel
 * 3) - 8) source region to be searched and resolution
 * 9) uncertainty (in percent from 0. to 1.) of the input group vel model
 * 10) safe half width (in second) of the surface wave signal (use a larger width if not sure)
 */
int main(int argc, char *argv[])
{
   if (argc!=11) {
      std::cout<<"Usage: "<<argv[0]<<" [Input SAC envelope list (fsac cper (optional grvel) )] [Input_Map (lon lat value)] [lonmin] [lonmax] [lonstep] [latmin] [latmax] [latstep] [grvel uncertainty (0.-1.)] [half_width of SW signal (in sec)]"<<std::endl;
      return -1;
   }

   /* Input map name */
   char mapname[100]; sprintf(mapname, "%s", argv[2]);
   /* Output map boundaries and grids */
   const float lonmin = atof(argv[3]), lonmax = atof(argv[4]), lonstep = atof(argv[5]);
   const float latmin = atof(argv[6]), latmax = atof(argv[7]), latstep = atof(argv[8]);
   const float safefactor = 1+atof(argv[9]), safehfwidth = atof(argv[10]);
   Point<float> LL(lonmin, latmin), UR(lonmax, latmax);
   Point<float> loccur, P1, P2;
   /* Output data matrix built with vector */
   const int nlon = (int)ceil((lonmax-lonmin)/lonstep)+1, nlat = (int)ceil((latmax-latmin)/latstep)+1; 
   std::vector< std::vector<float> > Probability_add(nlon, std::vector<float>(nlat, 0.));
   std::vector< std::vector<float> > Probability_add_all(nlon, std::vector<float>(nlat, 0.));
   std::vector< std::vector<float> > Probability_mul(nlon, std::vector<float>(nlat, 1.));
   std::vector< std::vector<float> > Probability_mul_all(nlon, std::vector<float>(nlat, 1.));
   std::vector< std::vector<int> > Counter_all(nlon, std::vector<int>(nlat, 0)), Counter(nlon, std::vector<int>(nlat, 0));
   //Probability.resize(nlon);
   //for(int i=0; i<nlon; i++) Probability[i].resize(nlat);
   

   /* read in and analyse each signal envelope */
   char fsac[150], cstaOld[10] = {"nan\0"};
   float cper, grv, lamda;
   float precamp, prectime;
   SACREC *sigcur = NULL;
   Path pathcur;
   /* read in the file list */
   std::fstream fsaclst;
   fsaclst.open(argv[1], std::ios::in);
   std::vector<std::string> filebuff;
   char buff[300];
   int ifile;
   for(ifile=0; fsaclst.getline(buff, 300); ifile++) filebuff.push_back(buff);
   fsaclst.close();
   int nfile = ifile;
   /* work on each file */
   bool SkipCurrentCsta = false;
   for(ifile=0; ifile<nfile; ifile++) {
      /* read in sac-file-name, center-period, and group-speed */
      grv = -1.;
      if( sscanf(filebuff.at(ifile).c_str(), "%s %f %f", fsac, &cper, &grv) < 2 ) {
			std::cerr<<"Warning(main): format error in file "<<argv[1]<<std::endl;
			continue;
      }
      std::cout<<"Searching precursors from file "<<fsac<<std::endl;
      /* load sac header and signal */
      sigcur = new SACREC(fsac, cper, grv);
      if( ! sigcur->load() ) { 
			delete sigcur; sigcur = NULL;
			std::cerr<<"Warning(main): sac file not exist OR bad data in sac file "<<fsac<<std::endl;
			continue;
      }
      P1 = sigcur->P1();
      P2 = sigcur->P2();
      /* predict grv of the current path if not provided */
      if( grv < 0. ) {
			new(&pathcur) Path(mapname, P1, P2);
			grv = pathcur.PathAverage(cper);
			if( grv < 0. ) { 
				std::cerr<<"Warning(main): meaningless group vel computed ("<<grv<<") for fsac "<<fsac<<std::endl;
				delete sigcur; sigcur = NULL; continue; 
			}
			sigcur->setGrv(grv);
			pathcur.~Path();
      }

      lamda = cper * grv * 0.5;
      /* checkpoint: print out a map for each center station */
      if( strcmp(cstaOld, "nan")==0 ) { sprintf(cstaOld, "%s", sigcur->Sta1()); }
      else if( strcmp(sigcur->Sta1(), cstaOld)==0 ) { 
			if(SkipCurrentCsta) { 
				std::cout<<"   data read from the hard drive already. skip..."<<std::endl; 
				delete sigcur; sigcur = NULL; continue; 
			}
      }
      else if( SkipCurrentCsta ) { SkipCurrentCsta = false; sprintf(cstaOld, "%s", sigcur->Sta1()); }// std::cerr<<"skip writing last csta"<<std::endl; }
      else {
	/* write out results for the last center station */
	 char outname[300];
	 std::fstream fout;
	 int i, j;
	 /* output Probabilit add and mul */
	 sprintf( outname, "Precursor_Probability_%s", cstaOld );
	 fout.open(outname, std::ios::out);
	 for( i=0,loccur=LL; loccur.Lon()<lonmax; loccur.move(lonstep, 0.),i++ ) for( j=0,loccur.SetLat(latmin); loccur.Lat()<latmax; loccur.move(0., latstep),j++ ) {
	    fout<<loccur.Lon()<<" "<<loccur.Lat()<<" "<<Probability_add[i][j]<<" "<<Probability_mul[i][j]<<" "<<Counter[i][j]<<std::endl;
	    Probability_add_all[i][j] += Probability_add[i][j]; Probability_add[i][j] = 0.;
	    Probability_mul_all[i][j] *= Probability_mul[i][j]; Probability_mul[i][j] = 1.;
	    Counter_all[i][j] += Counter[i][j]; Counter[i][j] = 0;
	 }
	 fout.close();
         //std::cerr<<"Checkpoint at ifile = "<<ifile<<" with (lon lat Padd Pmul) followed:"<<std::endl;
	 sprintf(cstaOld, "%s", sigcur->Sta1());
      }

      /* Check if a result file for the current center station exists. Read in the data and skip processing if it does */
      {
	 char inname[300];
	 sprintf( inname, "Precursor_Probability_%s", sigcur->Sta1() );
	 std::ifstream fin;
	 fin.open(inname);
	 /* if file exist. read in data and check for consistency */
	 if( fin ) {
	    int i, j;
	    float lontmp, lattmp;
	    char buff[300];
	    bool filematch = true;
	    std::vector< std::vector<float> > Padd(nlon, std::vector<float>(nlat, 0.)), Pmul(nlon, std::vector<float>(nlat, 1.));
	    std::vector< std::vector<int> > Ctr(nlon, std::vector<int>(nlat, 0));
            for( i=0,loccur=LL; loccur.Lon()<lonmax; loccur.move(lonstep, 0.),i++ ) {
	       if( ! filematch ) break;
	       for( j=0,loccur.SetLat(latmin); loccur.Lat()<latmax; loccur.move(0., latstep),j++ ) {
		  fin.getline(buff, 300);
		   sscanf(buff, "%f %f %f %f %d", &lontmp, &lattmp, &(Padd[i][j]), &(Pmul[i][j]), &(Ctr[i][j]) );
		   if( fabs(lontmp-loccur.Lon())>1.e-5 || fabs(lattmp-loccur.Lat())>1.e-5 ) { filematch = false; break; }
               }
	    }
	    fin.close();
	    /* if data looks good, merge them into the main matrices and skip all related paths */
	    if( filematch ) {
		for( i=0,loccur=LL; loccur.Lon()<lonmax; loccur.move(lonstep, 0.),i++ ) for( j=0,loccur.SetLat(latmin); loccur.Lat()<latmax; loccur.move(0., latstep),j++ )
		   { Probability_add_all[i][j] += Padd[i][j]; Probability_mul_all[i][j] *= Pmul[i][j]; Counter_all[i][j]+=Ctr[i][j];}
		SkipCurrentCsta = true;
		delete sigcur; sigcur = NULL;
		continue;
	    }
	 }

      }


      /* get location of the 2 stations */
      int i, j, idone = 0;
      float ftmp;
      #pragma omp parallel for private(j, loccur, pathcur)
      for( i=0; i<nlon; i++ ) {
	 loccur.SetLon(lonmin + i*lonstep);
	 for( j=0,loccur.SetLat(latmin); loccur.Lat()<latmax; loccur.move(0., latstep),j++ ) {
	    // predict travel time to sta1
	    new(&pathcur) Path(mapname, loccur, P1);
	    double t1 = pathcur.Dist() / pathcur.PathAverage(lamda); 
	    pathcur.~Path();
	    // predict travel time to sta2
	    new(&pathcur) Path(mapname, loccur, P2);
	    double t2 = pathcur.Dist() / pathcur.PathAverage(lamda); 
	    pathcur.~Path();
	    // check t1 and t2 for errors
	    if( t1<0. || t2<0. ) continue;
	
	    /* amp of precursor signal at time t2 - t1. Skip if merged with surface wave signal
	     * to be safe, both the uncertainty in the predicted arrival time and the half width of 
	     * the surface wave signal will be considered */
	    if( fabs(t2-t1) * safefactor > sigcur->Dist() / grv - safehfwidth ) continue;
	    #pragma omp critical
	    {
	    ftmp = sigcur->amp(t2-t1)/sigcur->noise();
		/* debug for nan
		if( ftmp != ftmp ) {
		   std::cerr<<"nan encountered working on path "<<fsac<<"("<<P1<<" - "<<P2<<") at sourcce loc "<<loccur<<std::endl;
		   std::cerr<<"t1="<<t1<<"  t2="<<t2<<"  sigamp="<<sigcur->amp(t2-t1)<<"  signoise="<<sigcur->noise()<<std::endl;
		   exit(0);
		}
		*/
	    Probability_add[i][j] += ftmp;
	    Probability_mul[i][j] *= ftmp;
	    Counter[i][j]++;
	    }
	    //std::cerr<<loccur<<": P="<<Probability_add[i][j]<<" T(precursor)="<<t2-t1<<" T(surface wave)="<<sigcur->Dist()/grv<<std::endl;
	 }
	 std::cout.precision(3);
        /* debug for nan
	#pragma omp critical
	{
	 idone++;
	 std::cout<<"   "<<ftmp*idone<<"% done from thread "<<omp_get_thread_num()<<"..."<<std::endl;
	}
	*/
      }
      delete sigcur; sigcur = NULL;
   }
   /* merge and write out results from the last center station */
   char outname[300];
   std::fstream fout;
   int i, j;
   /* output Probabilit add and mul */
   if( !SkipCurrentCsta ) {
      sprintf( outname, "Precursor_Probability_%s", cstaOld );
      fout.open(outname, std::ios::out);
      for( i=0,loccur=LL; loccur.Lon()<lonmax; loccur.move(lonstep, 0.),i++ ) for( j=0,loccur.SetLat(latmin); loccur.Lat()<latmax; loccur.move(0., latstep),j++ ) {
         fout<<loccur.Lon()<<" "<<loccur.Lat()<<" "<<Probability_add[i][j]<<" "<<Probability_mul[i][j]<<" "<<Counter[i][j]<<std::endl;
         Probability_add_all[i][j] += Probability_add[i][j]; 
         Probability_mul_all[i][j] *= Probability_mul[i][j]; 
	 Counter_all[i][j] += Counter[i][j];
      }
      fout.close();
   }


	

   /* output Probability add and mul*/
   sprintf(outname, "Precursor_Probability_all");
   fout.open(outname, std::ios::out);
   for( i=0,loccur=LL; loccur.Lon()<lonmax; loccur.move(lonstep, 0.),i++ )
      for( j=0,loccur.SetLat(latmin); loccur.Lat()<latmax; loccur.move(0., latstep),j++ ) fout<<loccur.Lon()<<" "<<loccur.Lat()<<" "<<Probability_add_all[i][j]<<" "<<Probability_mul_all[i][j]<<" "<<Counter_all[i][j]<<std::endl;
   fout.close();

   /* output Probability_mul */
//   sprintf(outname, "%s_Precursor_Probability_Mul", argv[1]);
//   fout.open(outname, std::ios::out);
//   for( i=0,loccur=LL; loccur.Lon()<lonmax; loccur.move(lonstep, 0.),i++ )
//      for( j=0,loccur.SetLat(latmin); loccur.Lat()<latmax; loccur.move(0., latstep),j++ ) fout<<loccur.Lon()<<" "<<loccur.Lat()<<" "<<Probability_mul_all[i][j]<<std::endl;
//   fout.close();

 return 0;
}
