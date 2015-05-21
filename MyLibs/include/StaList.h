#ifndef STALIST_H
#define STALIST_H

#ifndef FuncName
#define FuncName __FUNCTION__
#endif

#include "Point.h"
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <stdexcept>


struct StaInfo : public Point<float> {
   std::string name;
	/*
   //float lon, lat;
   static constexpr float maxmisloc = 0.001; // allow ~0.2km mislocation
   static constexpr float maxmislocS = maxmisloc*maxmisloc;
   static constexpr float NaN = -12345.;
	*/

public:
   StaInfo()
		: Point<float>() {}

   StaInfo( const std::string& namein, const float lonin, const float latin )
      : Point<float>(lonin, latin), name(namein) {}

	StaInfo( const std::string& line, const bool storeline = false ) { 
		if( ! LoadLine( line, storeline ) )
			throw std::runtime_error(std::string("Error(")+FuncName+"): ill-formed line "+line);
	}

   virtual bool LoadLine( const std::string& line, bool storeline = false ) {
		std::stringstream ss(line);
		bool suc = (ss >> lon >> lat >> name);
		if( ! suc ) {
			std::stringstream ss(line);
			suc = (ss >> name >> lon >> lat);
		}
		if( ! suc ) return false;
		if( lon < 0. ) lon += 360.;
		if( storeline ) name = line;

		return true;
   }

   bool IsSameName( const StaInfo& s2 ) const { return (name == s2.name); }
	/*
   bool IsSameLocation( const StaInfo& s2 ) const {
      float dislon = lon-s2.lon, dislat = lat-s2.lat;
      return (dislon*dislon + dislat*dislat < maxmislocS);
   }
	*/

   friend bool operator== ( const StaInfo& s1, const StaInfo& s2 ) {
      return ( s1.name==s2.name && s1.lon==s2.lon && s1.lat==s2.lat );
   }
   friend bool operator!= ( const StaInfo& s1, const StaInfo& s2 ) { return (! (s1==s2)); }

   friend std::ostream& operator<< ( std::ostream& o, const StaInfo& si ) {
      o<<si.name<<"  "<<std::setprecision(7)<<si.lon<<" "<<std::setprecision(7)<<si.lat;
      return o;
   }

};

class StaList {
public:
	StaList() {}
	StaList( const std::string& fname, const bool storeline = false ) { Load(fname, storeline); }

	void SortByLon() {
		std::sort( _dataV.begin(), _dataV.end() );
		sorted = true;
	}

	bool SearchLoc( const float lon, const float lat, StaInfo& data_find ) {
		if( ! sorted ) SortByLon();
      // iterator to the first data in _dataV with a lon>=loc.lon-StaInfo::maxmisloc
		StaInfo loc( "", lon, lat );
      auto loclb = loc; loclb.lon -= StaInfo::maxmisloc;
      auto idatal = std::lower_bound( _dataV.begin(), _dataV.end(), loclb );
      // iterator to the first data in dataV with a lon>loc.lon+StaInfo::maxmisloc
      auto locub = loc; locub.lon += StaInfo::maxmisloc;
      auto idatau = std::upper_bound( idatal, _dataV.end(), locub );
      auto idata = idatal;
		//std::cerr<<lon<<" "<<lat<<" "<<idatal-_dataV.begin()<<" "<<idatau-_dataV.begin()<<"\n";
      for( ; idata<idatau; idata++ )
         if( loc.IsSameLocation(*idata) ) break;
      if( idata < idatau ) { // found
			data_find = *idata;
			return true;
      } else { // not found
			return false;
      }
	}

	const std::vector<StaInfo>& DataRef() { return _dataV; }

	size_t size() { return _dataV.size(); }

private:
	std::vector<StaInfo> _dataV;
	bool sorted = false;

	void Load( const std::string& fname, const bool storeline = false ) {
		_dataV.clear();
		/* read in data file */
		std::ifstream fin(fname);
		if( ! fin )
			throw std::runtime_error(std::string("Error(")+FuncName+"): Cannot read from file "+fname);
		for(std::string line; std::getline(fin, line); ){
			try {
				_dataV.push_back( StaInfo(line, storeline) );
			} catch ( const std::exception& e ) {
				std::cerr<<e.what()<<"\n";
			}
		}
		//std::cout<<"   "<<dataV.size()<<" data lines read in."<<std::endl;
	}

};


#endif
