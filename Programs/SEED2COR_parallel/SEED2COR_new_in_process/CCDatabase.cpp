#include "CCDatabase.h"
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>

/*-------------------------------------- CCDatabase ----------------------------------------*/
/* Constructor: read in parameters -> fill in seed file list -> fill in station list */
CCDatabase::CCDatabase( const char *fname ) {
   CCParams.Load( fname );
   seedlst.Load( CCParams.seedfname );
   stalst.Load( CCParams.stafname );
	chlst.Load( CCParams.chlst_info );
	FillDInfo();
}

/* fill dinfo with parameters in CCParams */
void CCDatabase::FillDInfo() {
	dinfo.rdsexe = CCParams.rdsexe;
	dinfo.sps = CCParams.sps;
	dinfo.perl = CCParams.perl;
	dinfo.perh = CCParams.perh;
	dinfo.t1 = CCParams.t1;
	dinfo.tlen = CCParams.tlen;
	dinfo_rdy = false;
}

/* pull out the next daily record from the database */
bool CCDatabase::NextRecTest() {
   std::cerr<<*(seedlst.GetRec())<<std::endl;
   if( seedlst.ReLocate( 2012, 1, 26 ) ) { std::cerr<<*(seedlst.GetRec())<<std::endl; }
   if( seedlst.NextRec() ) std::cerr<<*(seedlst.GetRec())<<std::endl;
   if( seedlst.ReLocate( 2011, 2, 26 ) ) { std::cerr<<*(seedlst.GetRec())<<std::endl; }
   if( seedlst.NotEnded() ) std::cerr<<*(seedlst.GetRec())<<std::endl;
   if( seedlst.ReLocate( 2012, 2, 26 ) ) { std::cerr<<*(seedlst.GetRec())<<std::endl; }
   if( seedlst.NotEnded() ) std::cerr<<*(seedlst.GetRec())<<std::endl;

   std::cerr<<*(stalst.GetRec())<<std::endl;
   if( stalst.ReLocate( "J23A" ) ) { std::cerr<<*(stalst.GetRec())<<std::endl; }
   if( stalst.ReLocate( "SAO" ) ) { std::cerr<<*(stalst.GetRec())<<std::endl; }
   if( stalst.NextRec() ) std::cerr<<*(stalst.GetRec())<<std::endl;
   if( stalst.ReLocate( "CMA" ) ) { std::cerr<<*(stalst.GetRec())<<std::endl; }
   if( stalst.NotEnded() ) std::cerr<<*(stalst.GetRec())<<std::endl;

   return true;
}

bool CCDatabase::NextRec() {
	if( chlst.IsEnded()  ||
		 stalst.IsEnded() ||
		 seedlst.IsEnded() ) return false;
	dinfo_rdy = false;
	if( chlst.NextRec() ) return true;
	chlst.Rewind();
   if( stalst.NextRec() ) return true;
   stalst.Rewind();
   if( seedlst.NextRec() ) return true;
   return false;
}

void CCDatabase::Rewind() {
	chlst.Rewind();
	stalst.Rewind();
	seedlst.Rewind();
}

bool CCDatabase::GetRec(DailyInfo& dinfoout) {
	if( chlst.IsEnded()  ||
		 stalst.IsEnded() ||
		 seedlst.IsEnded() ) return false;
	if( ! dinfo_rdy ) {
		dinfo.Update( *(seedlst.GetRec()), *(stalst.GetRec()), *(chlst.GetRec()) );
		dinfo_rdy = true;
	}
	dinfoout = dinfo;
	return true;
}


/*-------------------------------------- CCPARAM ----------------------------------------*/
static int isTermi(int deno) {
   while(deno%2==0) deno *= 0.5;
   while(deno%5==0) deno *= 0.2;
   return deno==1;
}
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
		std::cerr<<"   Warning: cannot access rdseed through "<<rdsexe<<std::endl;
		if( FindInPath( "rdseed", rdsexe) )
			std::cerr<<"            corrected to "<<rdsexe<<std::endl;
		else
			ERR = 1;
   }
   //evrexe
   evrexe = (*fveciter).substr( 0, (*fveciter).find_first_of(" \t") ); fveciter++;
   std::cout<<"evalresp excutable:\t"<<evrexe<<std::endl;
   if( ! evrexe.find("evalresp") ) std::cout<<"   Warning: Are you sure this is an evalresp excutable?"<<std::endl;
   if( access(evrexe.c_str(), R_OK)!=0 ) {
		std::cerr<<"   Warning: cannot access evalresp through "<<evrexe<<std::endl;
		if( FindInPath( "evalresp", evrexe) )
         std::cerr<<"            corrected to "<<evrexe<<std::endl;
      else
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
	chlst_info = *(fveciter++);
   std::cout<<"channel list info:\t"<<chlst_info.substr(0, chlst_info.find("#"))<<std::endl;
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



/*-------------------------------------- Channellist ----------------------------------------*/
/* load in channel list from input channel info */
void Channellist::Load( const std::string& chinfo ) {
	std::stringstream ss( chinfo );
	std::string channel;
	while( (ss >> channel) && (channel.at(0) != '#') ) {
		icurrent = find(list.begin(), list.end(), channel );
		if( icurrent != list.end() ) {
			std::cerr<<"Warning(Channellist::Load): "<<channel<<" already in the list. Will be ignored!"<<std::endl;
			continue;
		}
		list.push_back(channel);
	}
	icurrent = list.begin();

	std::cout<<"Channellist::Load: "<<list.size()<<" channels loaded ( ";
	for(unsigned int i=0; i<list.size(); i++) std::cout<<list.at(i)<<" "; 
	std::cout<<")"<<std::endl;
}


/*-------------------------------------- Seedlist ----------------------------------------*/
/* load in seed list from input file and sort it by date */
static bool SameDate(const SeedInfo& a, const SeedInfo& b) {
   return ( a.year==b.year && a.month==b.month && a.day==b.day  );
}
static bool CompareDate(const SeedInfo& a, const SeedInfo& b) {
   return ( ( a.year<b.year ) || ( a.year==b.year && (a.month<b.month||(a.month==b.month&&a.day<b.day)) ) );
}
void Seedlist::Load( const std::string& fname ) {
	//list = new std::vector<SeedInfo>;
	std::ifstream fseed(fname);
	if( !fseed ) {
		std::cerr<<"ERROR(Seedlist::Load): Cannot open file "<<fname<<std::endl;
		exit(-1);
	}
	std::string buff;
	SeedInfo SRtmp;
	for(;std::getline(fseed, buff);) {
		char stmp[buff.length()];
		if( (sscanf(buff.c_str(),"%s %d %d %d", stmp, &(SRtmp.year), &(SRtmp.month), &(SRtmp.day))) != 4 ) { 
			std::cerr<<"Warning(Seedlist::Load): format error in file "<<fname<<std::endl; 
			continue;
		}
		SRtmp.name = stmp;
		list.push_back(SRtmp);
		//std::cerr<<list.back()<<std::endl;
	}
	fseed.close();
   std::cout<<"Seedlist::Load: "<<list.size()<<" seeds loaded"<<std::endl;
   /* sort the list by date */
   std::stable_sort(list.begin(), list.end(), CompareDate);
   //for(int i=0; i<list.size(); i++) std::cerr<<list.at(i)<<std::endl;
   icurrent = list.begin();
}

/* Move icurrent to (or to after) the next match of the input SeedInfo 
   assuming a sorted list and relocating using binary search */
bool Seedlist::ReLocate( int year, int month, int day ) { 
   SeedInfo srkey("", year, month, day);
   icurrent = std::lower_bound(list.begin(), list.end(), srkey, CompareDate );
   if( icurrent>=list.end() || icurrent<list.begin() ) return false;
   if( !SameDate(*icurrent, srkey) ) return false;
   return true;
}


/*-------------------------------------- Stationlist ----------------------------------------*/
/* load in station list from input file */
struct StaFinder {
   StaInfo a;
   StaFinder(StaInfo b) : a(b) {}
   bool operator()(StaInfo b) { return a.name.compare(b.name)==0; }
};
void Stationlist::Load( const std::string& fname ) {
   //list = new std::vector<StaInfo>;
   std::ifstream fsta(fname);
   if( !fsta ) {
      std::cerr<<"ERROR(Stationlist::Load): Cannot open file "<<fname<<std::endl;
      exit(-1);
   }
   std::string buff;
   StaInfo SRtmp;
	for(;std::getline(fsta, buff);) {
		char stmp[buff.length()];
		if( (sscanf(buff.c_str(),"%s %f %f", stmp, &(SRtmp.lon), &(SRtmp.lat))) != 3 ) { 
			std::cerr<<"Warning(Stationlist::Load): format error in file "<<fname<<std::endl; 
			continue;
		}
		SRtmp.name = stmp;
		icurrent = find_if(list.begin(), list.end(), StaFinder(SRtmp) );
		if( icurrent != list.end() ) {
			if( *icurrent == SRtmp ) { 
				std::cerr<<"Warning(Stationlist::Load): "<<SRtmp<<" already in the list. Will be ignored!"<<std::endl;
				continue;
			}
			else {
				std::cerr<<"Error(Stationlist::Load): station name confliction detected: "<<*icurrent<<" - "<<SRtmp<<std::endl;
				exit(0);
			}
		}
		list.push_back(SRtmp);
		//std::cerr<<list.back().fname<<" "<<list.back().lon<<" "<<list.back().lat<<std::endl;
	}
	fsta.close();
	std::cout<<"Stationlist::Load: "<<list.size()<<" stations loaded"<<std::endl;
	icurrent = list.begin();
}

/* Move icurrent to the next match of the input StaInfo 
	icurrent=.end() if no such match is found */
bool Stationlist::ReLocate( const std::string& staname ) { 
	StaInfo srkey(staname.c_str(), 0., 0.);
   icurrent = find_if(list.begin(), list.end(), StaFinder(srkey) );
   if( icurrent>=list.end() || icurrent<list.begin() ) return false;
   return true;
}

