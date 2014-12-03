/* given a data file whose colum1=lon and colum2=lat, pick out the lines with locations that exist in the given location list */
/* data file: (lon lat ...)
	location list: (lon lat ...) */

#include <cstring>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>

struct StaInfo {
   std::string name;
   float lon, lat;
	static constexpr float maxmisloc = 0.01;
	static constexpr float maxmislocS = maxmisloc*maxmisloc;
	static constexpr float NaN = -12345.;

public:
   StaInfo( const std::string namein = "", float lonin = NaN, float latin = NaN )
      : name(namein), lon(lonin), lat(latin) {}

   bool IsSameSta( const StaInfo& s2 ) const { return (name == s2.name); }
	bool IsSameLocation( const StaInfo& s2 ) const {
		float dislon = lon-s2.lon, dislat = lat-s2.lat;
		return (dislon*dislon + dislat*dislat < maxmislocS);
	}

   friend bool operator== ( const StaInfo& s1, const StaInfo& s2 ) {
		return ( s1.name==s2.name && s1.lon==s2.lon && s1.lat==s2.lat ); 
	}
   friend bool operator!= ( const StaInfo& s1, const StaInfo& s2 ) { return (! (s1==s2)); }

   friend std::ostream& operator<< ( std::ostream& o, const StaInfo& si ) { 
      o<<si.name<<"  "<<std::setprecision(7)<<si.lon<<" "<<std::setprecision(7)<<si.lat; 
      return o;
   }
};

bool CompareLon ( const StaInfo& s1, const StaInfo& s2 ) { return (s1.lon<s2.lon); }


int main(int argc, char* argv[]) {
   if( argc!=4 && argc!=5 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [data_infile] [location_lst] [data_outname] [keep loc #cols? (Y/N)]"<<std::endl;
      exit(-1);
   }
	/* output same #colums with location_lst? */
	bool samecol = true;
	if( argc == 5 )
		if( strcmp(argv[4],"Y") == 0 ) {
			samecol=true;
		} else {
			samecol = false;
		}

   /* read in data file */
   std::ifstream fin(argv[1]);
   if( ! fin ) {
      std::cerr<<"Error(main): Cannot read from file "<<argv[1]<<std::endl;
      exit(-2);
   }
   std::vector<StaInfo> dataV;
   for(std::string line; std::getline(fin, line); ){
		StaInfo datacur;
		datacur.name = line;
		if( sscanf(line.c_str(), "%f %f", &(datacur.lon), &(datacur.lat)) != 2 ) {
			std::cerr<<"Warning(main): wrong format detected in "<<argv[1]<<std::endl;
			continue;
		}
		if( datacur.lon < 0. ) datacur.lon += 360.;
		dataV.push_back(datacur);
   }
   fin.close(); fin.clear();
   std::cout<<"   "<<dataV.size()<<" data lines read in."<<std::endl;
	if( dataV.size() == 0 ) exit(-3);

   /* read in location list */
   fin.open(argv[2]);
   if( ! fin ) {
      std::cerr<<"Error(main): Cannot read from file "<<argv[2]<<std::endl;
      exit(-2);
   }
   std::vector<StaInfo> locV;
   for(std::string line; std::getline(fin, line); ){
		StaInfo loccur;
		if( sscanf(line.c_str(), "%f %f", &(loccur.lon), &(loccur.lat)) != 2 ) {
			std::cerr<<"Warning(main): wrong format detected in "<<argv[2]<<std::endl;
			continue;
		}
		if( loccur.lon < 0. ) loccur.lon += 360.;
		locV.push_back(loccur);
   }
   fin.close();
   std::cout<<"   "<<locV.size()<<" locations read in."<<std::endl;
	if( locV.size() == 0 ) exit(-3);

   /* sort datas by longitude */
	std::sort(dataV.begin(), dataV.end(), CompareLon);

	/* for each location, search in the dataV for the same location */
	std::ofstream fout(argv[3]);
	if( ! fout ) {
      std::cerr<<"Error(main): Cannot write to file "<<argv[3]<<std::endl;
      exit(-2);
	}
	int nmatch = 0;
	for(const auto& loc : locV) {
		// iterator to the first data in dataV with a lon>=loc.lon-StaInfo::maxmisloc
		auto loclb = loc; loclb.lon -= StaInfo::maxmisloc;
		auto idatal = std::lower_bound( dataV.begin(), dataV.end(), loclb, CompareLon );
		// iterator to the first data in dataV with a lon>loc.lon+StaInfo::maxmisloc
		auto locub = loc; locub.lon += StaInfo::maxmisloc;
		auto idatau = std::upper_bound( idatal, dataV.end(), locub, CompareLon );
		auto idata = idatal;
		for( ; idata<idatau; idata++ )
			if( loc.IsSameLocation(*idata) ) break;
		if( idata < idatau ) { // found
			fout<<(*idata).name<<"\n"; nmatch++;
		} else { // not found
			if( samecol ) fout<<loc.lon<<" "<<loc.lat<<std::endl;
		}
	}
	fout.close();
   std::cout<<"   "<<nmatch<<" matches found."<<std::endl;

   return 0;
}
