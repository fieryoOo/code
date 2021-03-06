#include <SacRec.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <deque>
#include <string>
#include <algorithm>

struct StaInfo {
   std::string name;
   float lon, lat;

public:
   StaInfo( const char* namein, float lonin, float latin )
      : name(namein), lon(lonin), lat(latin) {}

   bool IsSameSta( const StaInfo& s2 ) { return (name == s2.name); }

   friend bool operator== ( const StaInfo& s1, const StaInfo& s2 ) { return ( s1.name==s2.name && s1.lon==s2.lon && s1.lat==s2.lat ); }
   friend bool operator!= ( const StaInfo& s1, const StaInfo& s2 ) { return (! (s1==s2)); }

   friend std::ostream& operator<< ( std::ostream& o, const StaInfo& si ) { 
      o<<si.name<<"  "<<std::setprecision(7)<<si.lon<<" "<<std::setprecision(7)<<si.lat; 
      return o;
   }
};

bool CompareName ( const StaInfo& s1, const StaInfo& s2 ) { return (s1.name<s2.name); }


int main(int argc, char* argv[]) {
   if( argc != 3 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [sac_lst] [out_name]"<<std::endl;
      exit(-1);
   }

   /* read in sac list */
   std::ifstream fin(argv[1]);
   if( ! fin ) {
      std::cerr<<"Error(main): Cannot read from file "<<argv[1]<<std::endl;
      exit(0);
   }
   std::vector<std::string> saclst;
   for(std::string line; std::getline(fin, line); ){
      std::stringstream sin(line);
      std::string sacname;
      sin >> sacname;
      saclst.push_back(sacname);
   }
   fin.close();
   std::cout<<"   "<<saclst.size()<<" sac names read in."<<std::endl;

   /* read station info from each sac header and store into the station list */
   std::deque<StaInfo> stalst;
   for(int isac=0; isac<saclst.size(); isac++) {
      // load sac header
      SacRec sacrec( saclst[isac].c_str() );
		try {
			sacrec.LoadHD();
		} catch (std::exception& e) {
			std::cerr<<"Warning(main): LoadHD failed -- "<<e.what()<<std::endl;
			continue;
		}
      // form StaInfo
		std::string staid = sacrec.ntname() + " " + sacrec.stname();
      SAC_HD& shd = sacrec.shd;
      StaInfo SIcur(staid.c_str(), shd.stlo, shd.stla);
      // search for this same station in stalst
      std::deque<StaInfo>::iterator ista;
      ista = std::lower_bound( stalst.begin(), stalst.end(), SIcur, CompareName );
      if( ista!=stalst.end() && SIcur.IsSameSta(*ista) ) { // found station with same name!
		if( SIcur != *ista ) { // but different lon/lat
			std::cerr<<"Warning(main): conflicted location found for station "<<SIcur.name<<" ( "<<SIcur<<" : "<<*ista<<" )"<<std::endl;
		}
		continue;
      }
      //stalst.push_back(SIcur);
      stalst.insert(ista, SIcur);
		std::cout<<"..."<<100.0*isac/saclst.size()<<"\% done...\n\x1b[A";
   }
   std::cout<<"   "<<stalst.size()<<" station infos loaded."<<std::endl;

   /* sort station by alphabetic order */
   //stalst.sort(CompareName);

   /* write to file */
   std::ofstream fout(argv[2]);
   if( ! fout ) {
      std::cerr<<"Error(main): Cannot write to file "<<argv[2]<<std::endl;
      exit(0);
   }
   for(std::deque<StaInfo>::iterator ista=stalst.begin(); ista!=stalst.end(); ista++) fout<<*ista<<std::endl;
   fout.close();

   return 0;
}
