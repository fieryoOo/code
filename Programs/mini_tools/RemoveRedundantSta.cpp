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
   if( argc!=3 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [sta_list_in] [sta_list_out]"<<std::endl;
      exit(-1);
   }

   /* read in original station list */
   std::ifstream fin(argv[1]);
   if( ! fin ) {
      std::cerr<<"Error(main): Cannot read from file "<<argv[1]<<std::endl;
      exit(-2);
   }
   std::vector<StaInfo> staV;
   for(std::string line; std::getline(fin, line); ){
		StaInfo datacur;
		char buff[100];
		if( sscanf(line.c_str(), "%s %f %f", &(buff), &(datacur.lon), &(datacur.lat)) != 3 ) {
			std::cerr<<"Warning(main): wrong format detected in "<<argv[1]<<std::endl;
			continue;
		}
		datacur.name = buff;
		if( datacur.lon < 0. ) datacur.lon += 360.;
		staV.push_back(datacur);
   }
   std::cout<<"   "<<staV.size()<<" stations read in."<<std::endl;
	if( staV.size() == 0 ) exit(-3);

	/* brute-force search for now */
	for( auto ista1=staV.begin(); ista1<staV.end(); ista1++ ) {
		const auto& sta1 = *ista1;
		for( auto ista2=ista1+1; ista2<staV.end(); ) {
			const auto& sta2 = *ista2;
			if( sta1.name == sta2.name ) {
				if( sta1.IsSameLocation( sta2 ) ) {
					std::cout<<"   station ("<<(*ista2)<<") erased!"<<std::endl;
					ista2 = staV.erase(ista2);
				} else {
					std::cerr<<"Warning(main): name confliction detected for "<<sta1<<"  -  "<<sta2<<std::endl;
					ista2++;
				}
			} else {
				ista2++;
			}
		}
	}
   std::cout<<"   "<<staV.size()<<" stations left."<<std::endl;

	/* output */
	std::ofstream fout(argv[2]);
	for( const auto& sta : staV ) {
		fout<<sta<<std::endl;
	}


   return 0;
}
