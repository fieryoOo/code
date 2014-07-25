#include <iostream>
#include <string>
#include <vector>
/* seed record wrapper consists of filename, year, month, and day */
struct SeedInfo {
   std::string name;
   int year, month, day;
   SeedInfo() : year(0.), month(0.), day(0.) {}
   SeedInfo(const char* inname, const int& yy, const int& mm, const int& dd) { if(inname)name.assign(inname); year=yy; month=mm; day=dd; }
   friend std::ostream& operator<< (std::ostream& o, const SeedInfo& sr) { o<<"( "<<sr.name<<" "<<sr.year<<" "<<sr.month<<" "<<sr.day<<" )"; return o; }
};


/* station wrapper consists of stationname, longitude, and latitude */
struct StaInfo{
   std::string name;
   float lon, lat;
   StaInfo() : lon(0.), lat(0.) {}
   StaInfo(const char* inname, const float& lonin, const float& latin) { if(inname)name.assign(inname); lon=lonin; lat=latin; }
   friend bool operator== (StaInfo& a, StaInfo& b) { return ( a.name.compare(b.name)==0 && a.lon==b.lon && a.lat==b.lat ); }
   friend std::ostream& operator<< (std::ostream& o, const StaInfo& sr) { o<<"( "<<sr.name<<" "<<sr.lon<<" "<<sr.lat<<" "<<" )"; return o; }
};



struct DailyInfoData {
   SeedInfo seed;
   StaInfo sta;

   std::string rdsexe;
   std::string chname;
   std::string osac_outname, fsac_outname;
   std::string rec_outname, resp_outname;
   float sps, perl, perh, t1, tlen;

   //const std::string MonthName[13] {
   const std::vector<std::string> MonthName {
      "INVALID",
      "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
      "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
   };

};

int main() {
	DailyInfoData did1;
	DailyInfoData did2(did1);
	did2 = did1;
	return 0;
}
