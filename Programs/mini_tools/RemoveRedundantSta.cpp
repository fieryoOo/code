#include <cstring>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>

struct StaInfo {
   std::string name;
   std::string line;
   float lon, lat;
	static constexpr float maxmisloc = 0.001;	// allow ~0.1km mislocation
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
		return ( s1.name==s2.name && s1.IsSameLocation( s2 ) ); 
	}
   friend bool operator!= ( const StaInfo& s1, const StaInfo& s2 ) { return (! (s1==s2)); }

   friend std::ostream& operator<< ( std::ostream& o, const StaInfo& si ) { 
      o<<si.name<<"  "<<std::setprecision(7)<<si.lon<<" "<<std::setprecision(7)<<si.lat; 
      return o;
   }
};

bool CompareLon ( const StaInfo& s1, const StaInfo& s2 ) { return (s1.lon<s2.lon); }


int main(int argc, char* argv[]) {
   if( argc!=4 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [sta_list_in] [sta_list_out] [remove confliction? (0=retain;1=remove;2=retain the first)]"<<std::endl;
      exit(-1);
   }

   /* read in original station list */
   std::ifstream fin(argv[1]);
   if( ! fin ) {
      std::cerr<<"Error(main): Cannot read from file "<<argv[1]<<std::endl;
      exit(-2);
   }
   std::vector<StaInfo> staV;
   for(StaInfo datacur; std::getline(fin, datacur.line); ){
		char buff[100];
		if( sscanf(datacur.line.c_str(), "%s %f %f", &(buff), &(datacur.lon), &(datacur.lat)) != 3 ) {
			std::cerr<<"Warning(main): wrong format detected in "<<argv[1]<<std::endl;
			continue;
		}
		datacur.name = buff;
		if( datacur.lon < 0. ) datacur.lon += 360.;
		staV.push_back(datacur);
   }
   std::cout<<"   "<<staV.size()<<" stations read in."<<std::endl;
	if( staV.size() == 0 ) exit(-3);

	/* confliction handling? */
	int hconf = atoi(argv[3]);
	bool remove1st = (hconf==1);
	bool remove2nd = (hconf!=0);

	/* brute-force search for now */
	for( auto ista1=staV.begin(); ista1<staV.end(); ) {
		const auto& sta1 = *ista1;
		bool cflctDetected = false;
		for( auto ista2=ista1+1; ista2<staV.end(); ) {
			const auto& sta2 = *ista2;
			if( sta1.name == sta2.name ) {	// same sta name
				// check for confliction
				bool rmsta2 = remove2nd;	// remove_conflictions or
				if( sta1.IsSameLocation( sta2 ) ) {	// same location
					rmsta2 = true;	// remove the second sta
					std::cout<<"redundant station "<<sta1<<"  -  "<<sta2<<std::endl;
				} else {
					cflctDetected = true;
					std::cerr<<"Warning(main): name confliction detected for "<<sta1<<"  -  "<<sta2<<std::endl;
				}
				// remove station 2
				if( rmsta2 ) {
               std::cout<<"   station 2 ("<<(*ista2)<<") erased!"<<std::endl;
               ista2 = staV.erase(ista2);
				} else {
					ista2++;
				}
			} else {
				ista2++;
			}
		}
		if( cflctDetected && remove1st ) {
			std::cout<<"   station 1 ("<<(*ista1)<<") erased!"<<std::endl;
			ista1 = staV.erase(ista1);
		} else {
			ista1++;
		}
	}
   std::cout<<"   "<<staV.size()<<" stations left."<<std::endl;

	/* output */
	std::ofstream fout(argv[2]);
	for( const auto& sta : staV ) {
		fout<<sta.line<<std::endl;
	}


   return 0;
}
