#include "DailyRec.h"
#include "SysTools.h"
#include <iostream>
#include <iomanip>
#include <fstream>
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
struct DailyRec::DRimpl {
   //-------------------------------------------------------------parameters----------------------------------------------------------//
   std::string rdsexe;			//rdseed excutable
   std::string evrexe;			//evalresp excutable
   std::string fseed;			//input seed file name
   std::string fosac;			//output original sac name
   std::string ffsac;			//output resped sac name
   std::string staname;			//station list (J23A -129.6825 +44.844)
   std::string chname;			//channel to be extracted
   std::string outdir;			//output directory
   int sps;				//target sampling rate
   float gapfrac;			//maximum allowed gap fraction in input sac record
   float t1;				//cutting begining time in sec
   float tlen;				//time-record length in sec
   float perl, perh;			//signal period band
   int tnorm_flag;			//temperal normalization method. (1 for onebit, 2 for running average, 3 for earthquake cutting)
   float Eperl, Eperh;			//estimated earthquake period band (no effect when tnorm_flag==1; set Eperl=-1 to turn off the filter)
   float timehlen;			//half len of time window for running average (temperal normalization) (no effect when tnorm_flag==1)
   float frechlen;			//half len of frec window for whitenning. recommended: 0.0002; 0: no norm; -1: use input-smoothing file.
   std::string fwname;			//input spectrum reshaping signal file (takes effect only when frechlen==-1)
   //---------------------------------------------------------------------------------------------------------------------------------//
   
   std::string fresp;			//extracted resp file name
   std::string tdir;			//temp dir for running rdseed and evalresp

   float lon, lat;
   int year, month, day;
   //int npts;
   //float t0, delta;

public:
   inline int isTermi(int deno) {
      while(deno%2==0) deno /= 2;
      while(deno%5==0) deno /= 5;
      return deno==1;
   }

   bool CheckExistence( SacRec& sacrec ) {
      //check for sac file
      //SacRec sacrec(fosac.c_str());
      /* load header */
      std::string sacname(outdir);
      sacname.append("/"); sacname.append(fosac);
      if( ! ( sacrec.LoadHD(sacname.c_str()) ) ) return false;
      //SAC_HD *shd = read_shd(fosac);
      //if( shd==NULL ) return 0;
      //npts = sacrec.shd.npts;
      //t0 = sacrec.AbsTime();
      //delta = sacrec.shd.delta;
      int nlist;
      //check for RESP file
      char respname[150];//, dir[150];
      //sprintf(dir, "%s/%s", sdb->mo[imonth].name, sdb->ev[ne].name);
      sprintf(respname, "RESP.*.%s.*.%s", staname.c_str(), chname.c_str());
      //char *list_result = List(outdir.c_str(), respname, 0, &nlist);
      std::vector<std::string> list_result;
      if( ! List(outdir.c_str(), respname, 0, list_result) ) return false;
      //if(list_result==NULL) return false;
      if(nlist>1) std::cerr << "Warning(CheckExistence): more than one RESP files found for station "<<staname.c_str()<<" on "<<year<<"/"<<month<<"/"<<day<<std::endl; //reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: more than one RESP files found ***");
      //std::istringstream list( list_result );
      //list >> fresp;
      //free(list_result);
      fresp = list_result.at(0);
      return true;
   }

   bool Seed2Sac( const char *tdir, std::vector<std::string> &filelist ) {
      unsigned timeseed = std::chrono::system_clock::now().time_since_epoch().count();
      std::default_random_engine generator (timeseed);
      std::uniform_real_distribution<float> distribution(0., 1.);
      auto rand = std::bind ( distribution, generator );

/*
      char tdir[100];
      int ithread = 1;
      sprintf(tdir, "./Working_Thread_%d", ithread);
*/

      char fname[100], str[300];
      sprintf(fname, "%s/from_seed", tdir);
      FILE *ff = fopen(fname,"w");
      fprintf (ff, "%s <<END\n", rdsexe.c_str());
      fprintf (ff, "%s\n", fseed.c_str());
      fprintf(ff,"\n");				/* out file */
      fprintf(ff,"\n");				/* volume */
      fprintf(ff,"d\n");			/* option */
      fprintf(ff,"\n");				/* summary file */
      fprintf(ff,"%s\n", staname.c_str() );	/* station list */
      fprintf(ff,"%s\n", chname.c_str() );	/* channel list */
      fprintf(ff,"\n");				/* network list */
      fprintf(ff,"\n");				/* Loc Ids */
      fprintf(ff,"1\n");			/* out format */
      fprintf(ff,"N\n");			// new version!!!!!!!!!!
      fprintf(ff,"N\n");			/* Output poles & zeroes */
      fprintf(ff,"0\n");			/* Check Reversal */
      fprintf(ff,"\n");				/* Select Data Type */
      fprintf(ff,"\n");				/* Start Time */
      fprintf(ff,"\n");				/* End Time */
      fprintf(ff,"\n");				/* Sample Buffer Length  */
      fprintf(ff,"Y\n");			/* Extract Responses */
      fprintf(ff,"quit\n");
      fprintf(ff,"END\n");
      fclose(ff);

      //extract SAC&RESP and mv into thread working directory
      //pthread_mutex_lock(&rdslock); //lock for rdseed and shell operations
      sprintf(str,"sh %s >& /dev/null", fname);
      std::vector<std::string> list_result;
      #pragma omp critical
      {
      system(str);

      /*---------- mv response file -----------*/   
      sprintf(str, "RESP.*.%s.*.%s", staname.c_str(), chname.c_str());
      int nlist = 0;
      //list RESP files in the current depth
      //char *list_result = List(".", str, 0, &nlist);
      if( List(".", str, 0, list_result) ) {
	 //pthread_mutex_unlock(&rdslock);
	 //return false;
      //}
      /*
      if(list_result==NULL) {
	 //pthread_mutex_unlock(&rdslock);
	 return false;
      }
      */
         char list_name[150];
         int offset, curp = 0;
/*
      while( (sscanf(&list[curp], "%s%n", list_name, &offset)) == 1 ) {
	 if(curp==0) {
	    resp_fname = new char[150];
	    sprintf(resp_fname,"%s/%s", outdir.c_str(), list_name);
	    Move(list_name, resp_fname);
	 }
	 else fRemove(list_name);
	 curp += offset;
      }
*/
      /*--------------take the 1st resp name----------------*/
      /*
      std::stringstream list( list_result );
      std::string stmp;
      list >> stmp;
      fresp = outdir;
      fresp.append("/"); fresp.append(stmp);
      */
         fresp = outdir + "/" + list_result.at(0);
         Move(list_result.at(0).c_str(), fresp.c_str());
      /*--------------and remove the rest----------------*/
      //while( list >> stmp ) fRemove(stmp.c_str());
      //free(list_result);
         for(int i=1; i<list_result.size(); i++) fRemove(list_result.at(i).c_str());

      /*---------mv sac files and produce saclst---------*/
      //list SAC files
         sprintf(str, "*%s*%s*SAC", staname.c_str(), chname.c_str());
      }
      } // critical section
   
      //pthread_mutex_unlock(&rdslock); //unlock
      if( list_result.size() == 0 ) return false;
      return wMove(".", str, tdir, filelist);
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
DailyRec::DailyRec( const char *fname )
   : pimpl(new DRimpl()) { 
   if( fname ) Load(fname); 
}

DailyRec::DailyRec( const DailyRec& DRin )
   : pimpl(new DRimpl(*DRin.pimpl)) {}

DailyRec& DailyRec::operator= ( const DailyRec& DRin ) {
   pimpl.reset( new DRimpl(*DRin.pimpl) );
}

DailyRec::~DailyRec() {}//{ dRemove(pimpl->tdir.c_str()); }



/*---------------------------------------------------- Set/Load in parameters ----------------------------------------------------*/
int DailyRec::Set( const char *input ) {
   std::istringstream buff(input);
   std::string stmp;
   if( ! (buff>>stmp) ) return -2;
   bool succeed;
   if( stmp == "rdsexe" ) succeed = buff >> pimpl->rdsexe;
   else if( stmp == "evrexe" ) succeed = buff >> pimpl->evrexe;
   else if( stmp == "fseed" ) succeed = buff >> pimpl->fseed;
   else if( stmp == "fosac" ) succeed = buff >> pimpl->fosac;
   else if( stmp == "ffsac" ) succeed = buff >> pimpl->ffsac;
   else if( stmp == "staname" ) succeed = buff >> pimpl->staname;
   else if( stmp == "chname" ) succeed = buff >> pimpl->chname;
   else if( stmp == "outdir" ) succeed = buff >> pimpl->outdir;
   else if( stmp == "sps" ) succeed = buff >> pimpl->sps;
   else if( stmp == "gapfrac" ) succeed = buff >> pimpl->gapfrac;
   else if( stmp == "t1" ) succeed = buff >> pimpl->t1;
   else if( stmp == "tlen" ) succeed = buff >> pimpl->tlen;
   else if( stmp == "perl" ) succeed = buff >> pimpl->perl;
   else if( stmp == "perh" ) succeed = buff >> pimpl->perh;
   else if( stmp == "tnorm_flag" ) succeed = buff >> pimpl->tnorm_flag;
   else if( stmp == "Eperl" ) succeed = buff >> pimpl->Eperl;
   else if( stmp == "Eperh" ) succeed = buff >> pimpl->Eperh;
   else if( stmp == "timehlen" ) succeed = buff >> pimpl->timehlen;
   else if( stmp == "frechlen" ) succeed = buff >> pimpl->frechlen;
   else if( stmp == "fwname" ) succeed = buff >> pimpl->fwname;
   else return -1;

   if( succeed ) return 1;
   return 0;
}

bool DailyRec::Load( const char *fname ) {
   /* open param file */
   std::ifstream fparam(fname);
   if( !fparam ) {
      //std::cerr<<"ERROR(DailyRec::Load): Cannot open file "<<fname<<std::endl;
      return false;
   }
   /* read line by line for parameters */
   //std::cout<<"### Loading parameters from input file... ###"<<std::endl;
   int nparam = 0;
   std::string stmp;
   while( std::getline(fparam, stmp) ) {
      int retval = Set( stmp.c_str() );
      if( retval == -2 ) continue; // empty input
      else if( retval == -1 ){}// std::cerr<<"Warning(DailyRec::Load): Unknown parameter name: "<<stmp<<std::endl;
      else if( retval == 0 ){}// std::cerr<<"Warning(DailyRec::Load): Empty parameter field for "<<stmp<<std::endl;
      else nparam++;
   }
   fparam.close();
   std::cout<<"### "<<nparam<<" parameters loaded from param file "<<fname<<". ###"<<std::endl;
   return true;
}

/*---------------------------------------------------- Parameter Pre-Checks ----------------------------------------------------*/
bool DailyRec::CheckPreExtract() {
   /* Check each parameters. Return true only if all params look good */
   //std::cout<<"### Checking parameters for extracting sac... ###"<<std::endl;
   // rdsexe
   if( access( pimpl->rdsexe.c_str(), R_OK ) == -1 ) {
      std::cerr<<"Error(DailyRec::CheckPreExtract): Cannot access rdseed excutable at "<<pimpl->rdsexe.c_str()<<std::endl;
      return false;
   }
   // evrexe
   if( access( pimpl->evrexe.c_str(), R_OK ) == -1 ) {
      std::cerr<<"Error(DailyRec::CheckPreExtract): Cannot access evalresp excutable at "<<pimpl->evrexe.c_str()<<std::endl;
      return false;
   }
   // fseed
   if( access( pimpl->fseed.c_str(), R_OK ) == -1 ) {
      std::cerr<<"Error(DailyRec::CheckPreExtract): Cannot access input seed file at "<<pimpl->fseed.c_str()<<std::endl;
      return false;
   }
   // fosac
   if( pimpl->fosac.empty() ) {
      std::cerr<<"Error(DailyRec::CheckPreExtract): output osac name is empty"<<std::endl;
      return false;
   }
   // staname
   if( pimpl->staname.empty() ) {
      std::cerr<<"Error(DailyRec::CheckPreExtract): station name is empty"<<std::endl;
      return false;
   }
   // chname
   if( pimpl->chname.empty() ) {
      std::cerr<<"Error(DailyRec::CheckPreExtract): channel name is empty"<<std::endl;
      return false;
   }
   // outdir
   if( pimpl->outdir.empty() ) {
      std::cerr<<"Error(DailyRec::CheckPreExtract): output directory name is empty"<<std::endl;
      return false;
   }
   // sps
   if( ! ( pimpl->sps>0 && pimpl->isTermi(pimpl->sps) ) ) {
      std::cerr<<"Error(DailyRec::CheckPreExtract): 1/"<<pimpl->sps<<" isn't a terminating decimal"<<std::endl;
      return false;
   }
   // gapfrac
   if( pimpl->gapfrac<=0. || pimpl->gapfrac>1. ) {
      std::cerr<<"Error(DailyRec::CheckPreExtract): gapfrac("<<pimpl->gapfrac<<") is out of the range (0,1] "<<std::endl;
      return false;
   }

   std::cout<<"### Done(DailyRec::CheckPreExtract): all parameters checked. ###"<<std::endl;
   return true;
}

bool DailyRec::CheckPreRmRESP() {
   // ffsac
   if( pimpl->ffsac.empty() ) {
      std::cerr<<"Error(DailyRec::CheckPreRmRESP): output fsac name is empty"<<std::endl;
      return false;
   }
   // outdir
   if( pimpl->outdir.empty() ) {
      std::cerr<<"Error(DailyRec::CheckPreRmRESP): output directory name is empty"<<std::endl;
      return false;
   }
   // t1
   if( fabs(pimpl->t1) > 86400 ) std::cerr<<"Warning(DailyRec::CheckPreRmRESP): "<<pimpl->t1<<"sec exceeds one day."<<std::endl;
   // tlen
   float ftmp = pimpl->t1+pimpl->tlen;
   if( ftmp>86400 || ftmp<0 ) std::cerr<<"Warning(DailyRec::CheckPreRmRESP): ending time '"<<ftmp<<"sec' out of range."<<std::endl;
   // perl perh
   if( pimpl->perl<=0 || pimpl->perh<=0 ) {
      std::cerr<<"Error(DailyRec::CheckPreRmRESP): period band can not fall below 0."<<std::endl;
      return false;
   }
   if(pimpl->perl<2./pimpl->sps) {
      std::cerr<<"Error(DailyRec::CheckPreRmRESP): "<<pimpl->perl<<"sec is out of the lower limit at sps "<<pimpl->sps<<std::endl;
      return false;
   }
   if(pimpl->perl<5./pimpl->sps) std::cerr<<"Warning(DailyRec::CheckPreRmRESP): signal at "<<pimpl->perl<<"sec might be lost at a sampling rate of "<<pimpl->sps<<std::endl;

   std::cout<<"### Done(DailyRec::CheckPreRmRESP): all parameters checked. ###"<<std::endl;
   return true;
}

bool DailyRec::CheckPreTSNorm() {
   // tnorm_flag
   if( pimpl->tnorm_flag<0 || pimpl->tnorm_flag>4 ) {
      std::cerr<<"Error(DailyRec::CheckPreTSNorm): Unknow temperal norm method. Integer between 0 - 4 is expected"<<std::endl;
      return false;
   }
   // Eperl Eperh
   if(pimpl->tnorm_flag!=1) {
      if( pimpl->Eperl == -1. ) {}//std::cout<<"Eqk filter:\t\toff"<<std::endl;
      else {
         //cout<<"Eqk filter:\t\t"<<Eperl<<" - "<<Eperh<<"sec"<<endl;
         if( pimpl->Eperl<=0. || pimpl->Eperh<=0. ) {
            std::cerr<<"Error(DailyRec::CheckPreTSNorm): period band can not fall below 0."<<std::endl;
            return false;
         }
      }
   }
   // timehlen
   if(pimpl->tnorm_flag!=1) {
      //std::cout<<"t-len for run-avg:\t"<<timehlen<<std::endl;
      if(pimpl->timehlen<0.) {
         std::cerr<<"Error(DailyRec::CheckPreTSNorm): expecting a positive number for timehlen"<<std::endl;
         return false;
      }
   }
   // frechlen
   if(pimpl->frechlen==-1) {}//std::cout<<"input smoothing file will be used";
   else if(pimpl->frechlen<0.) {
      std::cerr<<"Error(DailyRec::CheckPreTSNorm): expecting a non-negative number for frechlen"<<std::endl;
      return false;
   }
   // fwname
   if(pimpl->frechlen==-1) {
      //std::cout<<"spec reshaping file:\t"<<pimpl->fwname<<std::endl;
      if( access(pimpl->fwname.c_str(), R_OK) == -1 ) {
         std::cerr<<"Error(DailyRec::CheckPreTSNorm): Cannot access the spec reshaping file "<<pimpl->fwname<<std::endl;
         return false;
      }
   }

   std::cout<<"### Done(DailyRec::CheckPreTSNorm): all parameters checked. ###"<<std::endl;
   return true;
}


/*---------------------------------------------------- Extract osac from seed file ----------------------------------------------------*/
bool DailyRec::ExtractSac( int fskipesac, bool writeout ) {
   std::ostringstream reports( std::ios_base::app );
   //std::ostringstream reports;
   bool bltmp = extractSac( fskipesac, writeout, reports );
   if( bltmp ) std::cout<<reports.str();
   return bltmp;
}
bool DailyRec::extractSac( int fskipesac, bool writeout, std::ostringstream& reports ) {
   /* random number generator */
   unsigned timeseed = std::chrono::system_clock::now().time_since_epoch().count();
   std::default_random_engine generator (timeseed);
   std::uniform_real_distribution<float> distribution(0., 1.);
   auto rand = std::bind ( distribution, generator );

   /* set up working directory */
   //std::string tdir;
   pimpl->tdir = "Working_" + std::to_string(rand());
  #pragma omp critical
  {
   while( ! MKDir( pimpl->tdir.c_str() ) ) pimpl->tdir = "Working_" + std::to_string(rand());
  }
   setvbuf(stdout, NULL, _IOLBF, 0);

   /* check for file existence */
   reports << "   Extracting Sac for station "<<pimpl->staname<<" on "<<pimpl->year<<"/"<<pimpl->month<<"/"<<pimpl->day<<": ";
   if( fskipesac==2 || fskipesac==1 ) {
      bool flag = pimpl->CheckExistence( sacT );
      if(fskipesac==2) {
	 reports << " skipped! "<<std::endl;
	 dRemove(pimpl->tdir.c_str());
	 return true;
      }
      if(flag) {
	 reports << " skipped on existense! "<<std::endl;
	 dRemove(pimpl->tdir.c_str());
	 return true;
      }
   }

   /* extract sacs from the seed */
   std::vector<std::string> filelst;
   //char *filelst = pimpl->Seed2Sac( pimpl->tdir.c_str(), &nfile );
   //if(filelst==NULL) {
 
   if(! pimpl->Seed2Sac( pimpl->tdir.c_str(), filelst ) ) {
      reports << " sac record not found from "<<pimpl->fseed<<"! "<<std::endl;
      dRemove(pimpl->tdir.c_str());
      return false;
   }

   /* read sacfile names from the filelst,
    * resample and merge one at a time */
   /*
   char sacname[150];
   int offset, curp = 0;
   sscanf(&filelst[curp], "%s%n", sacname, &offset);
   curp += offset;
   */
   sacT.Load(filelst.at(0).c_str());
   sacT.Resample(pimpl->sps);
   fRemove(filelst.at(0).c_str());
   bool merged = false;
   //while( (sscanf(&filelst[curp], "%s%n", sacname, &offset)) == 1 ) {
      //curp += offset;
   for(int i=1; i<filelst.size(); i++) {
      SacRec sacnew(filelst.at(i).c_str());
      sacnew.Load();
      sacnew.Resample(pimpl->sps);
      sacT.merge(sacnew);
      fRemove(filelst.at(i).c_str());
      merged = true;
   }
   //delete [] filelst;

   /* arrange the signal and check for percentage of poles if merged */
   if( merged ) {
      int Nholes;
      if( writeout ) {
	 std::string recname = pimpl->fosac + "_rec1";
	 Nholes = sacT.arrange(recname.c_str());
      }
      else Nholes = sacT.arrange();
      if( (float)Nholes/(float)sacT.shd.npts > pimpl->gapfrac ) {
         //pimpl->npts = -1;
	 reports << " too many holes detected! "<<std::endl;
	 dRemove(pimpl->tdir.c_str());
	 fRemove(pimpl->fresp.c_str());
         return false;
      }
   }

   /* output the merged sac file */
   if( writeout ) {
      std::string outname = pimpl->outdir + "/" + pimpl->fosac;
      sacT.Write(outname.c_str());
   }

   //dRemove(pimpl->tdir.c_str());
   reports << " done. "<<std::endl;
   dRemove(pimpl->tdir.c_str());
   return true;

/*
   if( MakeRecord(iev, ist, ithread) ) {
      if( nst%20 == 0 ) reports << "\n   "; //reports[ithread].tail += sprintf(reports[ithread].tail, "\n   ");
      //reports[ithread].tail += sprintf(reports[ithread].tail, "%s ", sdb->st[ist].name);
      reports << staname << " ";
      nst++;
   }
*/
}


/*---------------------------------------------------- Remove RESP ----------------------------------------------------*/

bool DailyRec::RmRESP(bool writeout) {
   if( pimpl->fresp.empty() ) return false;
   sacT.RmRESP(pimpl->evrexe.c_str(), pimpl->fresp.c_str(), pimpl->perl, pimpl->perh);
   fRemove(pimpl->fresp.c_str());
   // write out
   if( writeout ) {
      std::string outname = pimpl->outdir + "/" + pimpl->ffsac;
      sacT.Write(outname.c_str());
   }
   return true;
}


bool DailyRec::ZoomToEvent( std::string& ename, float evlon, float evlat, float tb, float te, bool writeout ) {
   if( te <= tb ) return false;
   
   // event time info
   SacRec sacevent;
   SAC_HD& eshd = sacevent.shd;
   eshd.nzyear = std::stoi( ename.substr(0, 4) );
   float month = std::stoi( ename.substr(4, 2) );
   float day = std::stoi( ename.substr(6, 2) );
   eshd.nzjday = pimpl->Jday( eshd.nzyear, month, day );
   eshd.nzhour = std::stoi( ename.substr(8, 2) );
   eshd.nzmin = std::stoi( ename.substr(10, 2) );
   eshd.nzsec = std::stoi( ename.substr(12, 2) );
   eshd.nzmsec = 0.;
   eshd.npts = 1; // so that AbsTime() won't fail

   // difference in origin time
   SAC_HD& shd = sacT.shd;
   float shift = sacevent.AbsTime() - sacT.AbsTime();
   // assign event time to sacrec and update shd.b
   shd.nzyear = eshd.nzyear; shd.nzjday = eshd.nzjday;
   shd.nzhour = eshd.nzhour; shd.nzmin = eshd.nzmin; shd.nzsec = eshd.nzsec;
   shd.b -= shift;

   // cut to tb - te
   if( ! sacT.cut( tb, te ) ) return false;

   // assign event location
   sprintf(shd.kevnm, "%s", ename.c_str());
   shd.evlo = evlon; shd.evla = evlat;

   // write to disk
   if( writeout ) {
      std::string outname = pimpl->outdir + "/" + pimpl->ffsac;
      sacT.Write(outname.c_str());
   }

   return true;
}
