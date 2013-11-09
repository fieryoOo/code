#include "PathAverage.h"
#include "SACREC.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <new>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_num_threads() { return 1;}
#endif

int main(int argc, char *argv[])
{
   if (argc!=9) {
      std::cout<<"Usage: "<<argv[0]<<" [Input SAC envelope list (fsac cper (optional grvel) )] [Input_Map (lon lat value)] [lonmin] [lonmax] [lonstep] [latmin] [latmax] [latstep]"<<std::endl;
      return -1;
   }

   /* Input map name */
   char mapname[100]; sprintf(mapname, "%s", argv[2]);
   /* Output map boundaries and grids */
   const float lonmin = atof(argv[3]), lonmax = atof(argv[4]), lonstep = atof(argv[5]);
   const float latmin = atof(argv[6]), latmax = atof(argv[7]), latstep = atof(argv[8]);
   Point<float> LL(lonmin, latmin), UR(lonmax, latmax);
   Point<float> loccur, P1, P2, P1old;
   /* Output data matrix built with vector */
   const int nlon = (int)ceil((lonmax-lonmin)/lonstep)+1, nlat = (int)ceil((latmax-latmin)/latstep)+1; 
   std::vector< std::vector<float> > Probability_add(nlon, std::vector<float>(nlat, 0.));
   std::vector< std::vector<float> > Probability_add_all(nlon, std::vector<float>(nlat, 0.));
   std::vector< std::vector<float> > Probability_mul(nlon, std::vector<float>(nlat, 1.));
   std::vector< std::vector<float> > Probability_mul_all(nlon, std::vector<float>(nlat, 1.));
   //Probability.resize(nlon);
   //for(int i=0; i<nlon; i++) Probability[i].resize(nlat);
   

   /* read in and analyse each signal envelope */
   int ifile;
   char buff[300], fsac[150];
   float cper, grv, lamda;
   float precamp, prectime;
   SACREC *sigcur = NULL;
   Path pathcur;
   std::fstream fsaclst;
   fsaclst.open(argv[1], std::ios::in);
   for(ifile=0; fsaclst.getline(buff, 300); ifile++) {
      /* read in sac-file-name, center-period, and group-speed */
      grv = -1.;
      if( sscanf(buff, "%s %f %f", fsac, &cper, &grv) < 2 ) {
	 std::cerr<<"Warning(main): format error in file "<<argv[1]<<std::endl;
	 continue;
      }
      std::cout<<"Searching precursors from file "<<fsac<<std::endl;
      /* load sac header and signal */
      sigcur = new SACREC(fsac, cper, grv);
      if( ! sigcur->load() ) { 
	 delete sigcur; sigcur = NULL;
	 std::cerr<<"Warning(main): bad data in sac file "<<fsac<<std::endl;
	 continue;
      }
      P1 = sigcur->P1();
      P2 = sigcur->P2();
      /* predict grv of the current path if not provided */
      if( grv < 0. ) {
	 new(&pathcur) Path(mapname, P1, P2);
	 grv = pathcur.PathAverage(cper);
	 if( grv < 0. ) continue;
	 sigcur->setGrv(grv);
	 pathcur.~Path();
      }

      lamda = cper * grv * 0.5;

      /* checkpoint: print out a map for each center station */
      if( ifile>0 && P1 != P1old ) {
	 char outname[300];
	 std::fstream fout;
	 int i, j;
	 /* output Probabilit add and mul */
	 sprintf( outname, "Precursor_Probability_%s", sigcur->Sta1() );
	 fout.open(outname, std::ios::out);
	 for( i=0,loccur=LL; loccur.Lon()<lonmax; loccur.move(lonstep, 0.),i++ ) for( j=0,loccur.SetLat(latmin); loccur.Lat()<latmax; loccur.move(0., latstep),j++ ) {
	    fout<<loccur.Lon()<<" "<<loccur.Lat()<<" "<<Probability_add[i][j]<<" "<<Probability_mul[i][j]<<std::endl;
	    Probability_add_all[i][j] += Probability_add[i][j]; Probability_add[i][j] = 0.;
	    Probability_mul_all[i][j] *= Probability_mul[i][j]; Probability_mul[i][j] = 1.;
	 }
	 fout.close();
         //std::cerr<<"Checkpoint at ifile = "<<ifile<<" with (lon lat Padd Pmul) followed:"<<std::endl;
	 P1old = P1;
      }

      /* get location of 2 the stations */
      /* search for largest precursoring signal */
      //if( ! sigcur->Precursor(&precamp, &prectime) ) continue;
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
	
	    // amp of precursor signal at time t2 - t1. Skip if merged with surface wave signal
	    // to be safe, a 5 percent uncertainty is assumed in the predicted arrival time
	    if( fabs(t2-t1) * 1.05 > sigcur->Dist() / grv - cper ) continue;
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
   fsaclst.close();

   char outname[300];
   std::fstream fout;
   int i, j;
   /* output Probability add and mul*/
   sprintf(outname, "Precursor_Probability_all");
   fout.open(outname, std::ios::out);
   for( i=0,loccur=LL; loccur.Lon()<lonmax; loccur.move(lonstep, 0.),i++ )
      for( j=0,loccur.SetLat(latmin); loccur.Lat()<latmax; loccur.move(0., latstep),j++ ) fout<<loccur.Lon()<<" "<<loccur.Lat()<<" "<<Probability_add_all[i][j]<<" "<<Probability_mul_all[i][j]<<std::endl;
   fout.close();

   /* output Probability_mul */
//   sprintf(outname, "%s_Precursor_Probability_Mul", argv[1]);
//   fout.open(outname, std::ios::out);
//   for( i=0,loccur=LL; loccur.Lon()<lonmax; loccur.move(lonstep, 0.),i++ )
//      for( j=0,loccur.SetLat(latmin); loccur.Lat()<latmax; loccur.move(0., latstep),j++ ) fout<<loccur.Lon()<<" "<<loccur.Lat()<<" "<<Probability_mul_all[i][j]<<std::endl;
//   fout.close();

 return 0;
}
