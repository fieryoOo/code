#ifndef CCDATABASE_H
#define CCDATABASE_H

#include <cstdio>
#include <vector>
#include <cstring>
#include <iostream>

/* ------------------------------ CCPARAM ------------------------------ */
/* CC parameter wrapper */
struct CCPARAM {
   //SAC_DB *sdb;
   int imonth;
   //-------------------------------------------------------------parameters----------------------------------------------------------//
   char rdsexe[200];		//rdseed excutable
   char evrexe[200];		//evalresp excutable
   char stafname[200];		//station list (J23A -129.6825 +44.844 optional_flag) (flag controls which sta-pairs to be ccd)
                        	// (skip cc between sta-pairs that are 1> both flagged 0 or 2> both not 0 but flagged with different numbers)
   char seedfname[200];		//SEED file list (down_US_ALL/OBS_BHZ_2012.FEB.29.203762.seed 2012 2 29)
   std::vector<std::string> chlst;//channel name list
   int sps;			//target sampling rate
   float gapfrac;		//maximum allowed gap fraction in input sac record
   float t1;			//cutting begining time in sec
   float tlen;			//time-record length in sec
   float perl, perh;		//signal period band
   int tnorm_flag;		//temperal normalization method. (1 for onebit, 2 for running average, 3 for earthquake cutting)
   float Eperl, Eperh;		//estimated earthquake period band (no effect when tnorm_flag==1; set Eperl=-1 to turn off the filter)
   float timehlen;		//half len of time window for running average (temperal normalization) (no effect when tnorm_flag==1)
   float frechlen;		//half len of frec window for whitenning. recommended: 0.0002; 0: no norm; -1: use input-smoothing file.
   char fwname[200];		//input spectrum reshaping signal file (takes effect only when frechlen==-1)
   int ftlen;			//turn on/off (1/0) cross-correlation-time-length correction for amplitude
   int fprcs;			//turn on/off (1/0) precursor signal checking
   float memomax;		//maximum memory fraction to be used. (set according to number of threads)
   int lagtime;			//cross-correlation signal half length in sec
   int mintlen;			//allowed minimum time length for cross-correlation (takes effect only when ftlen = 1)
   int fdelosac;		//delete original sac files after removing response when set to 1
   int fdelamph;		//delete am&ph files after cross-correlation when set to 1
   int fskipesac;		//skip ExtractSac() when set to 2, skip upon existence of target file when set to 1
   int fskipresp;		//skip RmRESP() when set to 2, skip upon existence of target file when set to 1
   int fskipamph;		//skip TempSpecNorm() when set to 2, skip upon existence of target file when set to 1
   int fskipcrco;		//skip CrossCorr() if target file exists when set to 1
   int CorOutflag;      	//controls the output of cross-corellation results: 0=monthly, 1=daily, 2=both
   //---------------------------------------------------------------------------------------------------------------------------------//
   CCPARAM() {}
   CCPARAM( const char* fname ) { Load(fname); }
   void Load( const char* );
};


/* ------------------------------ Seedlist ------------------------------ */
/* seed record wrapper consists of filename, year, month, and day */
struct SeedRec {
   char fname[200];
   int year, month, day;
   SeedRec() { year=0; month=0; day=0; }
   SeedRec(const char* inname, const int& yy, const int& mm, const int& dd) { if(inname)sprintf(fname, "%s", inname); year=yy; month=mm; day=dd; }
   friend std::ostream& operator<< (std::ostream& o, const SeedRec& sr) { o<<"( "<<sr.fname<<" "<<sr.year<<" "<<sr.month<<" "<<sr.day<<" )"; }
};
/* seed list with the 'NextRec' and the 'ReLocate' operation */
class Seedlist {
   std::vector<SeedRec> seedrec;
   std::vector<SeedRec>::iterator icurrent;
public:
   Seedlist() { icurrent = seedrec.end(); }
   /* call Load if fname is provided */
   Seedlist( const char* fname ) { Load(fname); }
   /* load station records from an input file */
   void Load( const char* );
   /* get to the next record and return true on success */
   bool NextRec( std::vector<SeedRec>::iterator& sriter) {
      if( icurrent<seedrec.begin() || icurrent>=seedrec.end() ) return false;
      sriter = icurrent++;
      return true; 
   }
   /* search for the first match of the input SeedRec and return true on success 
      icurrent will be moved to the match on succed and to the first rec after on failure */   
   bool ReLocate( int year, int month, int day );
   //~Seedlist() { if( seedrec ) delete seedrec; }
};


/* ------------------------------ Stationlist ------------------------------ */
/* station wrapper consists of stationname, longitude, and latitude */
struct StaRec{
   char fname[200];
   float lon, lat;
   StaRec() { lon=0.; lat=0.; }
   StaRec(const char* inname, const float& lonin, const float& latin) { if(inname)sprintf(fname, "%s", inname); lon=lonin; lat=latin; }
   friend bool operator== (StaRec& a, StaRec& b) { return ( strcmp(a.fname,b.fname)==0 && a.lon==b.lon && a.lat==b.lat ); }
   friend std::ostream& operator<< (std::ostream& o, const StaRec& sr) { o<<"( "<<sr.fname<<" "<<sr.lon<<" "<<sr.lat<<" "<<" )"; }
};
/* station list with the 'NextRec' and the 'ReLocate' operation */
class Stationlist {
   std::vector<StaRec> starec;
   std::vector<StaRec>::iterator icurrent;
public:
   Stationlist(){ icurrent = starec.begin(); }
   Stationlist( const char* fname ) { Load(fname); }
   /* load station records from an input file */
   void Load( const char* );
   /* get to the next record and return true on success */
   bool NextRec( std::vector<StaRec>::iterator& sriter) {
      if( icurrent<starec.begin() || icurrent>=starec.end() ) return false;
      sriter = icurrent++;
      return true; 
   }
   /* search for the first match of the input SeedRec and return true on success 
      icurrent=.end() if no such match is found */
   bool ReLocate( const char* staname );
};


/* ------------------------------ CCDatabase ------------------------------ */
class CCDatabase {
private:
   /* input parameters */
   CCPARAM CCParams;
   /* seed rec list */
   Seedlist seedlst;
   /* station list */
   Stationlist stalst;

public:
   /* constructor (read in parameters, seed list, and station list) */
   CCDatabase(char *inname);

   /* return a copy of the CC Parameters */
   const CCPARAM GetParams() { return CCParams; }

   /* Get the next daily record from the database. Assign file names and make directory if necessary */
   bool NextRec();
/*
   void InitialPthread();
   void FillMonths();
   void FillEvents();
*/

   //destructor (deallocate memories)
   ~CCDatabase(){};
//   void CleanupPthread();
};

#endif

