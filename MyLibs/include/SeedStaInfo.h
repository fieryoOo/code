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
	std::string osac_outname, fsac_outname;
	std::string rec_outname, resp_outname;
	float sps, perl, perh, t1, tlen;

	DailyInfo() {}
	DailyInfo( const SeedInfo sei, const StaInfo& sti )
		: SeedInfo(sei), StaInfo(sti) {}

	void Update( const SeedInfo sei, const StaInfo& sti ) {
		// copy from SeedInfo
		seedname = sei.seedname;
		year = sei.year; month = sei.month; day = sei.day;
		// copy from StaInfo
		staname = sti.staname;
		lon = sti.lon; lat = sti.lat;
		chname = "LHZ";
		osac_outname = "2012.SEP.1." + staname + "." + chname + ".SAC";
		fsac_outname = "ft_" + osac_outname;
		rec_outname = osac_outname + "_rec";
		resp_outname = "RESP.TA." + staname + ".." + chname;
	}

	friend std::ostream& operator<< (std::ostream& o, DailyInfo& di) { 
		o<<*(static_cast<SeedInfo*>(&di))<<"   "<<*(static_cast<StaInfo*>(&di)); 
		return o; 
	}

};


#endif
