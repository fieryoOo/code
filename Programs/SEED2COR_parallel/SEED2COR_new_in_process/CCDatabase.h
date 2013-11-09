#ifndef CCDATABASE_H
#define CCDATABASE_H

/* CC parameter wrapper */
struct CCPARAM {
   //SAC_DB *sdb;
   int imonth;
   //-------------------------------------------------------------parameters----------------------------------------------------------//
   char rdsexe[200];	//rdseed excutable
   char evrexe[200];	//evalresp excutable
   char stalst[200];	//station list (J23A -129.6825 +44.844 optional_flag) (flag controls which sta-pairs to be ccd)
                        // (skip cc between sta-pairs that are 1> both flagged 0 or 2> both not 0 but flagged with different numbers)
   char seedlst[200];	//SEED file list (down_US_ALL/OBS_BHZ_2012.FEB.29.203762.seed 2012 2 29)
   char ch[4];		//channel name
   int sps;		//target sampling rate
   float gapfrac;	//maximum allowed gap fraction in input sac record
   float t1;		//cutting begining time in sec
   float tlen;		//time-record length in sec
   float perl, perh;	//signal period band
   int tnorm_flag;	//temperal normalization method. (1 for onebit, 2 for running average, 3 for earthquake cutting)
   float Eperl, Eperh;	//estimated earthquake period band (no effect when tnorm_flag==1; set Eperl=-1 to turn off the filter)
   float timehlen;	//half len of time window for running average (temperal normalization) (no effect when tnorm_flag==1)
   float frechlen;	//half len of frec window for whitenning. recommended: 0.0002; 0: no norm; -1: use input-smoothing file.
   char fwname[200];	//input spectrum reshaping signal file (takes effect only when frechlen==-1)
   int ftlen;		//turn on/off (1/0) cross-correlation-time-length correction for amplitude
   int fprcs;		//turn on/off (1/0) precursor signal checking
   float memomax;	//maximum memory fraction to be used. (set according to number of threads)
   int lagtime;		//cross-correlation signal half length in sec
   int mintlen;		//allowed minimum time length for cross-correlation (takes effect only when ftlen = 1)
   int fdelosac;        //delete original sac files after removing response when set to 1
   int fdelamph;        //delete am&ph files after cross-correlation when set to 1
   int fskipesac;	//skip ExtractSac() when set to 2, skip upon existence of target file when set to 1
   int fskipresp;	//skip RmRESP() when set to 2, skip upon existence of target file when set to 1
   int fskipamph;	//skip TempSpecNorm() when set to 2, skip upon existence of target file when set to 1
   int fskipcrco;	//skip CrossCorr() if target file exists when set to 1
   int CorOutflag;      //controls the output of cross-corellation results: 0=monthly, 1=daily, 2=both
   //---------------------------------------------------------------------------------------------------------------------------------//
};

/* struct for seed list to be used by the CCDatabase class */
struct FSEED {
   int year, month, day;
   char name[200];
};


class CCDatabase {
private:
   /* input parameters */
   CCPARAM CCParams;
   /* seed records */
   FSEED fseed[];

public:
   /* constructor (read in parameters and station/seed list. Create DailyRecs for all stations all days) */
   CCDatabase(char *inname);

   /* return a copy of the CC Parameters */
   const CCPARAM GetParams() { return CCParams; }
/*
   int GetParameters(char *fname);
   int isTermi(int deno);
   void InitialPthread();
   void FillStations();
   void FillMonths();
   void FillEvents();

   // operations
   class DailyRec * PrepareNextRec(int *nmonth); // assign file names and make directory if necessary
*/
   //destructor (deallocate memories)
   ~CCDatabase(){};
//   void CleanupPthread();
};

#endif

