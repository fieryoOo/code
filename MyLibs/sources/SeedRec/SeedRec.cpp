#include "SysTools.h"
#include "SeedRec.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <unistd.h>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <functional>
//#include <omp.h>
#include "MyOMP.h"


/*---------------------------------------------------- implementation details ----------------------------------------------------*/
struct SeedRec::SRimpl {
   //-------------------------------------------------------------parameters----------------------------------------------------------//
   std::string rdsexe;			//rdseed excutable
   std::string fseed;			//input seed file name
   //---------------------------------------------------------------------------------------------------------------------------------//
 
   float lon, lat;
   int year, month, day;

	//static int num_obj;
	//int oid;

public:

	SRimpl( const std::string& fseedin, const std::string& rdsexein )
	: fseed(fseedin), rdsexe(rdsexein) {}

	bool FindInPath( const std::string fname, std::string& absname ) {
		char* pPath = getenv("PATH");
		if( ! pPath ) return false;
		std::stringstream sPath(pPath);
		for(std::string pathcur; std::getline(sPath, pathcur, ':'); ) {
			std::string testname = pathcur + '/' + fname;
			std::ifstream fin(testname);
			if( fin ) {
				absname = testname;
				return true;
			}
		}
		return false;
	}

   inline int isTermi(int deno) {
      while(deno%2==0) deno /= 2;
      while(deno%5==0) deno /= 5;
      return deno==1;
   }

   bool RunRdseed( const std::string& staname, const std::string& netname, const std::string& chname, const std::string& tdir,
						 std::string& fresp, std::vector<std::string> &filelist ) {

		// check parameters
		if( rdsexe.empty() || fseed.empty() || 
			 staname.empty() || chname.empty() ) return false;
		
		if( access( rdsexe.c_str(), R_OK ) == -1 ||
			 access( fseed.c_str(), R_OK ) == -1 ) return false;

      char fname[100], str[300];
      sprintf(fname, "%s/from_seed", tdir.c_str());
      FILE *ff = fopen(fname,"w");
      fprintf (ff, "%s <<- END\n", rdsexe.c_str());
      fprintf (ff, "%s\n", fseed.c_str());
      fprintf(ff,"\n");								/* out file */
      fprintf(ff,"\n");								/* volume */
      fprintf(ff,"d\n");							/* option */
      fprintf(ff,"\n");								/* summary file */
      fprintf(ff,"%s\n", staname.c_str() );	/* station list */
      fprintf(ff,"%s\n", chname.c_str() );	/* channel list */
		fprintf(ff,"\n");								/* network list */
		if( netname == "*" ) {						/* network list */
			fprintf(ff,"\n");
		} else {
			fprintf(ff,"%s\n", netname.c_str());
		}
      fprintf(ff,"\n");								/* Loc Ids */
      fprintf(ff,"1\n");							/* out format */
      fprintf(ff,"N\n");							/* new version!!!!!!!!!! */
      fprintf(ff,"N\n");							/* Output poles & zeroes */
      fprintf(ff,"0\n");							/* Check Reversal */
      fprintf(ff,"\n");								/* Select Data Type */
      fprintf(ff,"\n");								/* Start Time */
      fprintf(ff,"\n");								/* End Time */
      fprintf(ff,"\n");								/* Sample Buffer Length  */
      fprintf(ff,"Y\n");							/* Extract Responses */
      fprintf(ff,"quit\n");
      fprintf(ff,"END\n");
      fclose(ff);

      //extract SAC&RESP and mv into thread working directory
      sprintf(str,"sh %s >& /dev/null", fname);
      std::vector<std::string> list_result;
		#pragma omp critical(external)
      {
      system(str);

      /*---------- mv response file -----------*/   
      sprintf(str, "RESP.%s.%s.*.%s", netname.c_str(), staname.c_str(), chname.c_str());
      //list RESP files in the current depth
		fresp.clear();
      if( List(".", str, 0, list_result) ) {
         fresp = tdir + "/" + list_result.at(0);
         Move(list_result.at(0).c_str(), fresp.c_str());
			/*--------------and remove the rest----------------*/
         for(int i=1; i<list_result.size(); i++) fRemove(list_result.at(i).c_str());

      }
      } // critical section
   
      /*---------mv sac files and produce saclst---------*/
		//yyyy.ddd.hh.mm.ss.ffff.NN.SSSSS.LL.CCC.Q.SAC (rdseed naming convention; modify acordingly if the convention ever changes)
      sprintf(str, "*.*.*.*.*.*.*.%s.*.%s.*.SAC", staname.c_str(), chname.c_str());
      //if( list_result.size() == 0 ) return false;
      return wMove(".", str, tdir.c_str(), filelist);
   }

   int Jday ( int y, int m, int d ) {
      int i, jd = 0;
      for( i = 1; i < m; i++ ) {
	 if ( (i==1) || (i==3) || (i==5) || (i==7) || (i==8) || (i==10) ) jd += 31;
	 else if (i==2) {
	    if ( (y%400==0) || (y%100!=0&&y%4==0) ) jd += 29;
	    else jd += 28;
	 }
	 else jd += 30;
      }
      return jd + d;
   }

}; //pimpl


/*---------------------------------------------------- con/destructors & operations ----------------------------------------------------*/
SeedRec::SeedRec( const std::string& fname, std::ostream& reportin )
	: pimpl(new SRimpl(fname, "")), report(&reportin) {
	pimpl->FindInPath("rdseed", pimpl->rdsexe);
}
SeedRec::SeedRec( const std::string fname, const std::string rdsexein, std::ostream& reportin )
   : pimpl(new SRimpl(fname, rdsexein)), report(&reportin) {

	if( pimpl->rdsexe.empty() ) {
		pimpl->FindInPath("rdseed", pimpl->rdsexe);
	}

}

SeedRec::SeedRec( const SeedRec& SRin )
   : pimpl(new SRimpl(*SRin.pimpl)), report(SRin.report) {}

SeedRec::SeedRec( SeedRec&& SRin ) 
   : pimpl( std::move(SRin.pimpl) ), report(SRin.report) {}

SeedRec& SeedRec::operator= ( const SeedRec& SRin ) {
   pimpl.reset( new SRimpl(*SRin.pimpl) );
	report = SRin.report;
	return *this;
}

SeedRec& SeedRec::operator= ( SeedRec&& SRin ) {
   pimpl = std::move(SRin.pimpl);
	report = SRin.report;
   return *this;
}

SeedRec::~SeedRec() {}//{ dRemove(pimpl->tdir.c_str()); }



extern MEMO memo;
#include "MyOMP.h"
/*---------------------------------------------------- Extract osac from seed file ----------------------------------------------------*/
bool SeedRec::ExtractSac( const std::string& staname, const std::string& netname, const std::string& chname, const int sps,
								  const std::string& rec_outname, const std::string& resp_outname,
								  float& gapfrac, SacRec& sacout ) {
   /* random number generator */
   unsigned timeseed = std::chrono::system_clock::now().time_since_epoch().count();
   std::default_random_engine generator (timeseed);
   std::uniform_real_distribution<float> distribution(0., 1.);
   auto rand = std::bind ( distribution, generator );

   /* set up working directory */
   std::string tdir;
  #pragma omp critical (mkdir)
  {
   while( true ) {
		tdir = "Working_" + std::to_string(rand());
		if( MKDir( tdir.c_str() ) ) break;
	}
  }

   setvbuf(stdout, NULL, _IOLBF, 0);

   /* extract sacs from the seed */
	std::string fresp;
   std::vector<std::string> filelst;
   if( ! pimpl->RunRdseed( staname, netname, chname, tdir, fresp, filelst ) ) {
      //reports << " sac record not found from " << pimpl->fseed << "! " << std::endl;
      dRemove(tdir.c_str());
      return false;
   }

   /* read sacfile names from the filelst,
    * resample and merge one at a time */
   sacout.Load(filelst.at(0).c_str());
   sacout.Resample(sps);
   fRemove(filelst.at(0).c_str());
   bool merged = false;
   for(int i=1; i<filelst.size(); i++) {
      SacRec sacnew(filelst.at(i).c_str(), *report);
      sacnew.Load();
      sacnew.Resample(sps);
      sacout.merge(sacnew);
      fRemove(filelst.at(i).c_str());
      merged = true;
   }
   //delete [] filelst;

   /* arrange the signal and check for percentage of poles if merged */
	gapfrac = 0.;
   if( merged ) {
      int Nholes;
      if( ! rec_outname.empty() )
			Nholes = sacout.arrange(rec_outname.c_str());
      else 
			Nholes = sacout.arrange();
      gapfrac = (float)Nholes/(float)sacout.shd.npts;
   }

	// rename resp file
	if( ! resp_outname.empty() )
		Move( fresp.c_str(), resp_outname.c_str() );

   //reports << " done. "<<std::endl;
   dRemove(tdir.c_str());
   return true;

}


