#include "CCDatabase.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <algorithm>

/*-------------------------------------- CCDatabase ----------------------------------------*/
/* Constructor: read in parameters -> fill in seed file list -> fill in station list */
CCDatabase::CCDatabase( char *fname ) {
   CCParams.Load( fname );
   seedlst.Load( CCParams.seedfname.c_str() );
   stalst.Load( CCParams.stafname.c_str() );
}

/* pull out the next daily record from the database */
bool CCDatabase::NextRec() {
   std::vector<SeedRec>::iterator iseed;
   if( seedlst.NextRec( iseed ) ) std::cerr<<*iseed<<std::endl;
   if( seedlst.ReLocate( 2012, 1, 26 ) ) { seedlst.NextRec( iseed ); std::cerr<<*iseed<<std::endl; }
   if( seedlst.ReLocate( 2011, 2, 26 ) ) { seedlst.NextRec( iseed ); std::cerr<<*iseed<<std::endl; }

   std::vector<StaRec>::iterator ista;
   if( stalst.NextRec( ista) ) std::cerr<<*ista<<std::endl;
   if( stalst.ReLocate( "J23A" ) ) { stalst.NextRec( ista ); std::cerr<<*ista<<std::endl; }
   if( stalst.ReLocate( "SAO" ) ) { stalst.NextRec( ista ); std::cerr<<*ista<<std::endl; }
}


/*-------------------------------------- CCPARAM ----------------------------------------*/
static int isTermi(int deno) {
   while(deno%2==0) deno /= 2;
   while(deno%5==0) deno /= 5;
   return deno==1;
}
/* read in parameters for the CC Database from the inputfile */
void CCPARAM::Load( const char* fname ) {
   /* load param file input a vector */
   std::ifstream fparam(fname);
   if( ! fparam ) {
      std::cerr<<"Error(GetParam): Cannot open parameter file "<<fname<<std::endl;
      exit(0);
   }
   std::vector<std::string> filevec;
   for(std::string line; std::getline(fparam, line); ) filevec.push_back(line);
   fparam.close();
   if( filevec.size() < 29 ) { std::cerr<<"   Error(CCPARAM::Load): No enough param lines in file "<<fname<<std::endl; exit(0); }

   /* load/check one parameter at a time */
   int ERR=0, itmp;
   float ftmp;
   std::vector<std::string>::iterator fveciter = filevec.begin();
   std::cout<<"---------------------------Checking input parameters---------------------------"<<std::endl;
   //rdsexe
   rdsexe = (*fveciter).substr( 0, (*fveciter).find_first_of(" \t") ); fveciter++;
   std::cout<<"rdseed excutable:\t"<<rdsexe<<std::endl;
   if( ! rdsexe.find("rdseed") ) std::cout<<"   Warning: Are you sure this is an rdseed excutable?"<<std::endl;
   if( access(rdsexe.c_str(), R_OK)!=0 ) {
     std::cerr<<"   Error: cannot access rdseed through "<<rdsexe<<std::endl;
     ERR = 1;
   }
   //evrexe
   evrexe = (*fveciter).substr( 0, (*fveciter).find_first_of(" \t") ); fveciter++;
   std::cout<<"evalresp excutable:\t"<<evrexe<<std::endl;
   if( ! evrexe.find("evalresp") ) std::cout<<"   Warning: Are you sure this is an evalresp excutable?"<<std::endl;
   if( access(evrexe.c_str(), R_OK)!=0 ) {
     std::cerr<<"   Error: cannot access evalresp through "<<evrexe<<std::endl;
     ERR = 1;
   }
   //stafname
   stafname = (*fveciter).substr( 0, (*fveciter).find_first_of(" \t") ); fveciter++;
   std::cout<<"station list:\t\t"<<stafname<<std::endl;
   std::ifstream filetmp(stafname.c_str());
   if( ! filetmp ) {
      std::cerr<<"   Error: cannot access file "<<stafname<<std::endl;
      ERR = 1;
   }
   else {
      std::string buff;
      std::getline(filetmp, buff);
      char stmp[buff.length()];
      if( sscanf(buff.c_str(), "%s %f %f", stmp, &ftmp, &ftmp) != 3 ) {
         std::cerr<<"   Error: incorrect format in file "<<stafname<<"! Shoud be (sta lon lat)"<<std::endl;
         ERR = 1;
      }
      filetmp.close();
   }
   //seedfname
   seedfname = (*fveciter).substr( 0, (*fveciter).find_first_of(" \t") ); fveciter++;
   std::cout<<"seed list:\t\t"<<seedfname<<std::endl;
   filetmp.open(seedfname.c_str(), std::ifstream::in);
   if( ! filetmp ) {
      std::cerr<<"   Error: cannot access file "<<seedfname<<std::endl;
      ERR = 1;
   }
   else {
      std::string buff;
      std::getline(filetmp, buff);
      char stmp[buff.length()];
      if( sscanf(buff.c_str(), "%s %d %d %d", stmp, &itmp, &itmp, &itmp) != 4 ) {
         std::cerr<<"   Error: incorrect format in file "<<seedfname<<"! (path/seedfile yyyy mm dd) is expected"<<std::endl;
         ERR = 1;
      }
      filetmp.close();
   }
   //chlst
   std::stringstream ss( *(fveciter++) );
   std::string channel;
   while( (ss >> channel) && (channel.at(0) != '#') ) chlst.push_back(channel);
   std::cout<<"channel list:\t\t"<<chlst.size()<<" channels in the list ( ";
   for(int i=0; i<chlst.size(); i++) std::cout<<chlst.at(i)<<" "; std::cout<<")"<<std::endl;
   //sps
   sscanf((*(fveciter++)).c_str(), "%f", &ftmp);
   sps = static_cast<int>(ftmp);
   std::cout<<"target sampling rate:\t"<<ftmp<<std::endl;
   if( sps!=ftmp || sps <= 0 ) {
      std::cerr<<"   Error: a positive integer is expected!"<<std::endl;
      ERR = 1;
   }
   else if( !isTermi(sps) ) {
      std::cerr<<"   Error: 1/"<<sps<<" isn't a terminating decimal, will cause rounding error!"<<std::endl;
      ERR = 1;
   }
   //gapfrac
   sscanf((*(fveciter++)).c_str(), "%f", &gapfrac);
   std::cout<<"max gap fraction:\t"<<gapfrac*100<<"%"<<std::endl;
   if( gapfrac<0 || gapfrac>1 ) {
      std::cerr<<"   Error: a number between 0. - 1. is expected!"<<std::endl;
      ERR = 1;
   }
   //t1
   sscanf((*(fveciter++)).c_str(), "%f", &t1);
   std::cout<<"cutting begining:\t"<<t1<<"sec"<<std::endl;
   if( fabs(t1) > 86400 ) std::cout<<"   Warning: "<<t1<<"sec exceeds one day."<<std::endl;
   //tlen
   sscanf((*(fveciter++)).c_str(), "%f", &tlen);
   std::cout<<"time-rec length:\t"<<tlen<<"sec"<<std::endl;
   ftmp = t1+tlen;
   if( ftmp>86400 || ftmp<0 ) std::cout<<"   Warning: ending time '"<<ftmp<<"sec' out of range."<<std::endl;
   //perl perh
   sscanf((*(fveciter++)).c_str(), "%f", &perl);
   sscanf((*(fveciter++)).c_str(), "%f", &perh);
   std::cout<<"signal per band:\t"<<perl<<" - "<<perh<<"sec"<<std::endl;
   if( perl<0 || perh<0 ) {
      std::cerr<<"   Error: period band can not go below 0."<<std::endl;
      ERR = 1;
   }
   if(perl<2./sps) {
      std::cerr<<"   Error: "<<perl<<"sec is out of the lower limit at sps "<<sps<<std::endl;
      ERR = 1;
   }
   if(perl<5./sps) std::cout<<"   Warning: signal at "<<perl<<"sec might be lost at a sampling rate of "<<sps<<std::endl;
   //tnorm_flag
   sscanf((*(fveciter++)).c_str(), "%d", &tnorm_flag);
   std::cout<<"t-norm method:\t\t";
   if(tnorm_flag==0) std::cout<<"none"<<std::endl;
   else if(tnorm_flag==1) std::cout<<"One-bit"<<std::endl;
   else if(tnorm_flag==2) std::cout<<"Running average"<<std::endl;
   else if(tnorm_flag==3) std::cout<<"Earthquake cutting"<<std::endl;
   else {
      std::cerr<<std::endl<<"   Error: Unknow method. integer between 0 - 3 is expected"<<std::endl;
      ERR = 1;
   }
   //Eperl Eperh
   sscanf((*(fveciter++)).c_str(), "%f", &Eperl);
   sscanf((*(fveciter++)).c_str(), "%f", &Eperh);
   if(tnorm_flag!=1) {
      if( Eperl == -1. ) std::cout<<"Eqk filter:\t\toff"<<std::endl;
      else {
         std::cout<<"Eqk filter:\t\t"<<Eperl<<" - "<<Eperh<<"sec"<<std::endl;
         if( Eperl<0 || Eperh<0 ) {
            std::cerr<<"   Error: period band can not go below 0."<<std::endl;
            ERR = 1;
         }
      }
   }
   //timehlen
   sscanf((*(fveciter++)).c_str(), "%f", &timehlen);
   if(tnorm_flag!=1) {
      std::cout<<"t-len for run-avg:\t"<<timehlen<<std::endl;
      if(timehlen<0) {
         std::cerr<<"   Error: positive number is expected"<<std::endl;
         ERR = 1;
      }
   }
   //frechlen
   sscanf((*(fveciter++)).c_str(), "%f", &frechlen);
   std::cout<<"t-len for whitening:\t"<<frechlen<<"  ";
   if(frechlen==-1) std::cout<<"input smoothing file will be used";
   else if(frechlen<0) {
      std::cerr<<std::endl<<"   Error: non-negative number is expected";
      ERR = 1;
   }
   std::cout<<std::endl;
   //fwname
   fwname = (*fveciter).substr( 0, (*fveciter).find_first_of(" \t") ); fveciter++;
   if(frechlen==-1) {
      std::cout<<"spec reshaping file:\t"<<fwname<<std::endl;
      if( access(fwname.c_str(), R_OK)!=0 ) {
         std::cerr<<"   Error: cannot access file "<<fwname<<std::endl;
         ERR = 1;
      }  
   }
   //ftlen
   sscanf((*(fveciter++)).c_str(), "%d", &ftlen);
   if(ftlen<=0) std::cout<<"cor-time-len correction\toff"<<std::endl;
   else { std::cout<<"cor-time-len correction\ton"<<std::endl; /*if(tnorm_flag==3) ftlen = 2;*/ }
   //fprcs
   sscanf((*(fveciter++)).c_str(), "%d", &fprcs);
   if(fprcs<=0) std::cout<<"prcsr-signal checking\toff"<<std::endl;
   else std::cout<<"prcsr-signal checking\ton"<<std::endl;
   //memomax
   sscanf((*(fveciter++)).c_str(), "%f", &memomax);
   std::cout<<"memory fraction:\t"<<memomax*100<<"%"<<std::endl;
   if( memomax<0 || memomax>1 ) {
      std::cerr<<"   Error: a number between 0. - 1. is expected!"<<std::endl;
      ERR = 1;
   }
   //lagtime
   sscanf((*(fveciter++)).c_str(), "%d", &lagtime);
   std::cout<<"cor lag time:\t\t"<<lagtime<<"sec"<<std::endl;
   if(lagtime<0) {
      std::cerr<<"   Error: negative lagtime!"<<std::endl;
      ERR = 1;
   }
   else if(lagtime>524288) std::cout<<"   Warning: lag time exceeds the maximum"<<std::endl;
   //mintlen
   sscanf((*(fveciter++)).c_str(), "%d", &mintlen);
   std::cout<<"min time len:\t\t"<<mintlen<<"sec"<<std::endl;
   if( mintlen>86400 ) std::cout<<"   Warning: allowed minimum time length larger than a day."<<std::endl;
   //fdelosac
   sscanf((*(fveciter++)).c_str(), "%d", &fdelosac);
   std::cout<<"Delete orig sacs?\t";
   if(fdelosac==1) std::cout<<"Yes"<<std::endl;
   else std::cout<<"No"<<std::endl;
   //fdelamph
   sscanf((*(fveciter++)).c_str(), "%d", &fdelamph);
   std::cout<<"Delete am&ph files?\t";
   if(fdelamph==1) std::cout<<"Yes"<<std::endl;
   else std::cout<<"No"<<std::endl;
   //fskipesac
   sscanf((*(fveciter++)).c_str(), "%d", &fskipesac);
   std::cout<<"ExtractSac()\t\t";
   if(fskipesac==2) std::cout<<"skip"<<std::endl;
   else if(fskipesac==1) std::cout<<"skip if target file exists"<<std::endl;
   else std::cout<<"overwrite if target file exists"<<std::endl;
   //fskipresp
   sscanf((*(fveciter++)).c_str(), "%d", &fskipresp);
   std::cout<<"RmRESP()\t\t";
   if(fskipresp==2) std::cout<<"skip"<<std::endl;
   else if(fskipresp==1) std::cout<<"skip if target file exists"<<std::endl;
   else std::cout<<"overwrite if target file exists"<<std::endl;
   //fskipamph
   sscanf((*(fveciter++)).c_str(), "%d", &fskipamph);
   std::cout<<"TempSpecNorm()\t\t";
   if(fskipamph==2) std::cout<<"skip"<<std::endl;
   else if(fskipamph==1) std::cout<<"skip if target file exists"<<std::endl;
   else std::cout<<"overwrite if target file exists"<<std::endl;
   //fskipcrco
   sscanf((*(fveciter++)).c_str(), "%d", &fskipcrco);
   std::cout<<"CrossCorr()\t\t";
   if(fskipcrco==2) std::cout<<"skip"<<std::endl;
   else if(fskipcrco==1) std::cout<<"skip if target file exists"<<std::endl;
   else std::cout<<"overwrite if target file exists"<<std::endl;
   //CorOutflag
   sscanf((*(fveciter++)).c_str(), "%d", &CorOutflag);
   std::cout<<"CC Output: \t\t";
   switch( CorOutflag ) {
      case 0:
	 std::cout<<"monthly"<<std::endl;
	 break;
      case 1:
	 std::cout<<"daily"<<std::endl;
	 break;
      case 2:
	 std::cout<<"monthly + daily"<<std::endl;
	 break;
      default:
	 std::cerr<<std::endl<<"   Error: unknow outflag("<<CorOutflag<<")!"<<std::endl;
	 ERR = 1;
   }
   //check for EOF
   if( fveciter != filevec.end() ) std::cout<<"   Warning: End of file not reached!"<<std::endl;
   std::cout<<"-------------------------------Checking completed-------------------------------"<<std::endl;
   if(ERR==1) exit(0);
   /* prompt to continue */
   char cin;
   std::string line;
   std::cout<<"Continue?  ";
   std::getline(std::cin, line);
   sscanf(line.c_str(), "%c", &cin);
   if( cin!='Y' && cin!='y' ) exit(0);

}



/*-------------------------------------- Seedlist ----------------------------------------*/
/* load in seed list from input file and sort it by date */
static bool SameDate(const SeedRec& a, const SeedRec& b) {
   return ( a.year==b.year && a.month==b.month && a.day==b.day  );
}
static bool CompareDate(const SeedRec& a, const SeedRec& b) {
   return ( ( a.year<b.year ) || ( a.year==b.year && (a.month<b.month||(a.month==b.month&&a.day<b.day)) ) );
}
void Seedlist::Load( const char* fname ) {
   //seedrec = new std::vector<SeedRec>;
   std::ifstream fseed(fname);
   if( !fseed ) {
      std::cerr<<"ERROR(Seedlist::Load): Cannot open file "<<fname<<std::endl;
      exit(-1);
   }
   std::string buff;
   SeedRec SRtmp;
   for(;std::getline(fseed, buff);) {
      char stmp[buff.length()];
      if( (sscanf(buff.c_str(),"%s %d %d %d", stmp, &(SRtmp.year), &(SRtmp.month), &(SRtmp.day))) != 4 ) { 
	 std::cerr<<"Warning(Seedlist::Load): format error in file "<<fname<<std::endl; 
	 continue;
      }
      SRtmp.fname = stmp;
      seedrec.push_back(SRtmp);
      //std::cerr<<seedrec.back()<<std::endl;
   }
   fseed.close();
   std::cout<<"Seedlist::Load: "<<seedrec.size()<<" seeds loaded"<<std::endl;
   /* sort the list by date */
   std::stable_sort(seedrec.begin(), seedrec.end(), CompareDate);
   //for(int i=0; i<seedrec.size(); i++) std::cerr<<seedrec.at(i)<<std::endl;
   icurrent = seedrec.begin();
}

/* Move icurrent to (or to after) the next match of the input SeedRec 
   assuming a sorted list and relocating using binary search */
bool Seedlist::ReLocate( int year, int month, int day ) { 
   SeedRec srkey("", year, month, day);
   icurrent = std::lower_bound(seedrec.begin(), seedrec.end(), srkey, CompareDate );
   if( icurrent>=seedrec.end() || icurrent<seedrec.begin() ) return false;
   if( !SameDate(*icurrent, srkey) ) return false;
   return true;
}


/*-------------------------------------- Stationlist ----------------------------------------*/
/* load in station list from input file */
struct StaFinder {
   StaRec a;
   StaFinder(StaRec b) : a(b) {}
   bool operator()(StaRec b) { return a.fname.compare(b.fname)==0; }
};
void Stationlist::Load( const char* fname ) {
   //starec = new std::vector<StaRec>;
   std::ifstream fsta(fname);
   if( !fsta ) {
      std::cerr<<"ERROR(Stationlist::Load): Cannot open file "<<fname<<std::endl;
      exit(-1);
   }
   std::string buff;
   StaRec SRtmp;
   for(;std::getline(fsta, buff);) {
      char stmp[buff.length()];
      if( (sscanf(buff.c_str(),"%s %f %f", stmp, &(SRtmp.lon), &(SRtmp.lat))) != 3 ) { 
	 std::cerr<<"Warning(Stationlist::Load): format error in file "<<fname<<std::endl; 
	 continue;
      }
      SRtmp.fname = stmp;
      icurrent = find_if(starec.begin(), starec.end(), StaFinder(SRtmp) );
      if( icurrent != starec.end() ) {
	 if( *icurrent == SRtmp ) { 
	    std::cerr<<"Warning(Stationlist::Load): "<<SRtmp<<" already in the list. Will be ignored!"<<std::endl;
	    continue;
	 }
	 else {
	    std::cerr<<"Error(Stationlist::Load): station name confliction detected: "<<*icurrent<<" - "<<SRtmp<<std::endl;
	    exit(0);
	 }
      }
      starec.push_back(SRtmp);
      //std::cerr<<starec.back().fname<<" "<<starec.back().lon<<" "<<starec.back().lat<<std::endl;
   }
   fsta.close();
   std::cout<<"Stationlist::Load: "<<starec.size()<<" stations loaded"<<std::endl;
   icurrent = starec.begin();
}

/* Move icurrent to the next match of the input StaRec 
   icurrent=.end() if no such match is found */
bool Stationlist::ReLocate( const char* staname ) { 
   StaRec srkey(staname, 0., 0.);
   icurrent = find_if(starec.begin(), starec.end(), StaFinder(srkey) );
   if( icurrent>=starec.end() || icurrent<starec.begin() ) return false;
   return true;
}

