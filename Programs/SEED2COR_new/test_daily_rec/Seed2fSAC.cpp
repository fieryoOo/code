#include "SysTools.h"
#include "DailyRec.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <sstream>
//#include <omp.h>
#include "MyOMP.h"

struct EventInfo {
   std::string seed_name, event_name;
   float lon, lat;

   EventInfo( const char* evinfo ) {
      std::stringstream sin(evinfo);
      if( ! (sin >> seed_name >> event_name >> lon >> lat) ) {
	 std::cerr<<"Error(EventInfo::EventInfo): Incorrect input format ( "<<evinfo<<" )!"<<std::endl;
	 exit(0);
      }
   }

};

struct StaInfo {
   std::string name;
   float lon, lat;

public:
   StaInfo( const char* namein, float lonin, float latin )
      : name(namein), lon(lonin), lat(latin) {}

   StaInfo( const char* stainfo ) {
      std::stringstream sin(stainfo);
      if( ! (sin >> name >> lon >> lat) ) {
         std::cerr<<"Error(StaInfo::StaInfo): Incorrect input format ( "<<stainfo<<" )!"<<std::endl;
         exit(0);
      }
   }

   bool IsSameSta( StaInfo& s2 ) { return (name == s2.name); }

   friend bool operator== ( const StaInfo& s1, const StaInfo& s2 ) { return ( s1.name==s2.name && s1.lon==s2.lon && s1.lat==s2.lat ); }
   friend bool operator!= ( const StaInfo& s1, const StaInfo& s2 ) { return (! (s1==s2)); }

   friend std::ostream& operator<< ( std::ostream& o, const StaInfo& si ) {
      o<<si.name<<"  "<<si.lon<<" "<<si.lat;
      return o;
   }
};


int  main ( int argc, char *argv[] ) 
{
   if( argc != 4) {
      std::cerr<<"Usage: "<<argv[0]<<" [input parameter file] [seed_event_lst (seed_name event_name lon lat)] [station_lst]"<<std::endl;
      return -1;
   }

   /* read in seed_event list */
   std::ifstream fin(argv[2]);
   if( ! fin ) {
      std::cerr<<"Error(main): cannot read from file "<<argv[2]<<std::endl;
      exit(0);
   }
   std::vector<EventInfo> evlst;
   for(std::string line; std::getline(fin, line); ) {
      evlst.push_back( EventInfo( line.c_str() ) );
   }
   fin.close();
   std::cout<<"### "<<evlst.size()<<" seed_events read in. ###"<<std::endl;

   /* read in station list */
   fin.clear(); fin.open(argv[3]);
   if( ! fin ) {
      std::cerr<<"Error(main): cannot read from file "<<argv[3]<<std::endl;
      exit(0);
   }
   std::vector<StaInfo> stalst;
   for(std::string line; std::getline(fin, line); ) {
      stalst.push_back( StaInfo( line.c_str() ) );
   }
   fin.close();
   std::cout<<"### "<<stalst.size()<<" stations read in. ###"<<std::endl;


   /* store old SAC and RESP files if there's any */
   //bool oldfiles = false;
   MKDir("old_sac_files");
   //wMove(".", "*.SAC", "old_sac_files", 0, &nfmvd); if( nfmvd > 0 ) oldfiles = true;
   //wMove(".", "RESP.*", "old_sac_files", 0, &nfmvd); if( nfmvd > 0 ) oldfiles = true;
   std::vector<std::string> filelst;
   bool oldsac = wMove(".", "*.SAC", "old_sac_files", filelst);
   bool oldresp = wMove(".", "RESP.*", "old_sac_files", filelst);


   //int ithread; 
   /* create the DailyRec object, load in parameters from the input file */
   DailyRec dailyrec(argv[1]);

   /* channel list */
   std::vector<std::string> chlst = {"LHZ", "LHN", "LHE"};
   /* iterate through all events */
   for(int iev=0; iev<evlst.size(); iev++) {
      // copy params from dailyrec
      DailyRec DRtmp_ev(dailyrec);

      // set fseed and outdir for this event
      std::string stmp_ev = "fseed " + evlst[iev].seed_name;
      DRtmp_ev.Set(stmp_ev.c_str());
      stmp_ev = "outdir " + evlst[iev].event_name;
      DRtmp_ev.Set(stmp_ev.c_str());
      MKDir( evlst[iev].event_name.c_str() );

      // iterate through all stations
      #pragma omp parallel for schedule(dynamic, 1)
      for(int ista=0; ista<stalst.size(); ista++) {
	 // set staname
	 std::string stmp = "staname " + stalst[ista].name;
	 DailyRec DRtmp_sta( DRtmp_ev );
	 DRtmp_sta.Set(stmp.c_str());

	 for(int ich=0; ich<chlst.size(); ich++) {
	    DailyRec DRtmp( DRtmp_sta );
	    std::string& ename(evlst[iev].event_name);

	    // set fosac
	    stmp = "fosac " + ename + "." + stalst[ista].name + "." + chlst[ich] + ".sac";
	    DRtmp.Set(stmp.c_str());
	    // set ffsac
	    stmp = "ffsac " + ename + "." + stalst[ista].name + "." + chlst[ich] + ".sac";
	    DRtmp.Set(stmp.c_str());
	    // set chname
	    stmp = "chname " + chlst[ich];
	    DRtmp.Set(stmp.c_str());

	    //if( ! DRtmp.CheckPreExtract() ) exit(0);

	    // extract sac: donot skip | don't write yet
	    if( ! DRtmp.ExtractSac(0, false) ) continue;

	    // reResponse: remove RESP file | don't write yet
	    DRtmp.RmRESP(false);

	    // cut record by event: write to disk
	    DRtmp.ZoomToEvent( ename, evlst[iev].lon, evlst[iev].lat, 0, 5000, true );
	 }
      }
      std::cout<<"### event "<<evlst[iev].event_name<<" is done!"<<std::endl;
   }
   std::cout<<"### all done! ###"<<std::endl;
   /* fetch back old files and remove temporary dir */
   /*
   if( oldfiles ) {
      wMove("old_sac_files", "*.SAC", ".", 0, &nfmvd);
      wMove("old_sac_files", "RESP.*", ".", 0, &nfmvd);
   }
   */
   if( oldsac ) wMove("old_sac_files", "*.SAC", ".", filelst);
   if( oldresp ) wMove("old_sac_files", "RESP.*", ".", filelst);
   fRemove("old_sac_files");

   return 0;
}
