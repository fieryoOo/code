#ifndef SEEDSTAINFO_H
#define SEEDSTAINFO_H

#include <iostream>


/* seed record wrapper consists of filename, year, month, and day */
struct SeedInfo {
   std::string seedname;
   int year, month, day;
   SeedInfo() { year=0; month=0; day=0; }
   SeedInfo(const char* inname, const int& yy, const int& mm, const int& dd) { if(inname)seedname.assign(inname); year=yy; month=mm; day=dd; }
   friend std::ostream& operator<< (std::ostream& o, const SeedInfo& sr) { o<<"( "<<sr.seedname<<" "<<sr.year<<" "<<sr.month<<" "<<sr.day<<" )"; return o; }
};


/* station wrapper consists of stationname, longitude, and latitude */
struct StaInfo{
   std::string staname;
   float lon, lat;
   StaInfo() { lon=0.; lat=0.; }
   StaInfo(const char* inname, const float& lonin, const float& latin) { if(inname)staname.assign(inname); lon=lonin; lat=latin; }
   friend bool operator== (StaInfo& a, StaInfo& b) { return ( a.staname.compare(b.staname)==0 && a.lon==b.lon && a.lat==b.lat ); }
   friend std::ostream& operator<< (std::ostream& o, const StaInfo& sr) { o<<"( "<<sr.staname<<" "<<sr.lon<<" "<<sr.lat<<" "<<" )"; return o; }
};


/* Daily Info derived from StaInfo and SeedInfo */
struct DailyInfo : public SeedInfo, public StaInfo {
	std::string rdsexe;
	std::string chname;
	std::string rec_outname, resp_outname;
	float sps;

	DailyInfo( const SeedInfo sei, const StaInfo& sti, const std::string rdsexe )
		: SeedInfo(sei), StaInfo(sti), rdsexe(rdsexe)
		, chname("LHZ"), sps(1), rec_outname("TEST/temp_rec.txt"), resp_outname("TEST/RESP_temp") {}

	friend std::ostream& operator<< (std::ostream& o, DailyInfo& di) { 
		o<<*(static_cast<SeedInfo*>(&di))<<"   "<<*(static_cast<StaInfo*>(&di)); 
		return o; 
	}

};


#endif
