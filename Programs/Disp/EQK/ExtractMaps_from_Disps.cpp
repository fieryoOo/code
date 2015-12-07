/*
This code extracts information (distance, azimuth, traveltime, amplitude, and snr) from FTAN results of earthquakes.
( !!! Note that the input (dist_in_sac) and output distance could be different, whereas the input is expected to be the one 
used by the aFTAN (which is usually read from sac header) and the output is computed from the input source-station location !!! )
Input 1: station.lst (STANM LON LAT)
Input 2: input-file list (STA1 STA2 Disp_File SNR_File dist_in_sac)
*/
#include "DisAzi.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class PointC {
public:
   std::string name, net;
   float lon, lat;

public:
   /* con/destructors and operators */
   PointC() {}

   PointC( const char *info ) { Load(info); }

   PointC( const char *namein, float lonin, float latin )
      : name(namein), lon(lonin), lat(latin) {}

   ~PointC() {}

   //friend bool operator== ( PointC &a, PointC &b ) { return ( a.name==b.name && a.lon==b.lon && a.lat==b.lat );  }

   friend std::ostream& operator<< (std::ostream& o, const PointC& pt) {
      o<<pt.name<<" "<<pt.lon<<" "<<pt.lat;
      return o;
   }

   /* load point from string */
   void Load( const char *info ) {
      std::stringstream sin(info);
      if( ! (sin >> name >> lon >> lat) ) {
			sin.clear(); sin.str(info);
			sin >> net >> name >> lon >> lat;
		}
      if( lon < 0 ) lon += 360.;
   }
};

class PathRec {
   std::vector<PointC> *stalst;
   std::vector<PointC> *evelst;
   std::string fdisp, fsnr, famp;
   /* fill a point by searching the sta/eve list*/
   bool FillSrc( const char *name ) { return FillPointC(*evelst, name, src); }

   bool FillSta( const char *name ) { return FillPointC(*stalst, name, sta); }

   bool FillPointC( std::vector<PointC> &ptlst, const char *name, PointC &pt ) {
      int i;
      for(i=0; i<ptlst.size(); i++) {
			if( ptlst[i].name == name ) break;
		}
      if( i == ptlst.size() ) return false;
      pt = ptlst.at(i);
      return true;
   }

   void InserSort(float *arr, float *dat, int n) {
      int i,j;
      float ftmp0, ftmp1;
      for(i=1;i<n;i++) {
         ftmp0=arr[i]; ftmp1=dat[i];
         for(j=i;j>0 && ftmp0<arr[j-1]; j--) { arr[j]=arr[j-1]; dat[j]=dat[j-1]; }
         arr[j]=ftmp0; dat[j]=ftmp1;
      }
   }

   void InserSort2(float *arr, float *dat1, float *dat2, int n) {
      int i,j;
      float ftmp0, ftmp1, ftmp2;
      for(i=1;i<n;i++) {
         ftmp0=arr[i]; ftmp1=dat1[i]; ftmp2=dat2[i];
         for(j=i;j>0 && ftmp0<arr[j-1]; j--) { arr[j]=arr[j-1]; dat1[j]=dat1[j-1]; dat2[j]=dat2[j-1]; }
         arr[j]=ftmp0; dat1[j]=ftmp1; dat2[j]=ftmp2;
      }
   }

public:
   bool loaded, hasAmpfile = false;
   PointC src, sta;
   double distsac, dist, azi1, azi2;
   float snr, amp, grT, phT;

   /* con/destructors and operators */
   PathRec( std::vector<PointC> &stalstin, std::vector<PointC> &evelstin, const char *fdispin, const char *fsnrin, const char *srcname, const char *staname, double& distin )
      : stalst(&stalstin), evelst(&evelstin)
      , fdisp(fdispin), fsnr(fsnrin), distsac(distin) {
      Initialize();
      loaded = FillSrc( srcname ) && FillSta( staname );
   }

   const char* Format() {
      std::string format("evname(1) evlon(2) evlat(3)  staname(4) stalon(5) stalat(6)  dist(7) azi1(8) azi2(9) : snr(11) amp(12) grT(13) phT(14)");
      return format.c_str();
   }

   friend std::ostream& operator<< (std::ostream& o, const PathRec& path) {
      o<<path.src<<"  "<<path.sta<<"  "<<path.dist<<" "<<path.azi1<<" "<<path.azi2<<" : "<<path.snr<<" "<<path.amp<<" "<<path.grT<<" "<<path.phT;
      return o;
   }

   /* load in path from a info string */
   PathRec( std::vector<PointC> &stalstin, std::vector<PointC> &evelstin, const char *info )
      : stalst(&stalstin), evelst(&evelstin) {
      Initialize();
      std::stringstream sin(info);
      std::string staname, srcname;
      if( ! (sin >> srcname >> staname >> fdisp >> fsnr >> distsac) )
			throw std::runtime_error(std::string("Error(PathRec::PathRec): format error (")+info+").");
		hasAmpfile = (sin >> famp);
      loaded = FillSrc( srcname.c_str() ) && FillSta( staname.c_str() );
		//std::cerr<<"srcname = "<<srcname<<" staname = "<<staname<<" fdisp = "<<fdisp<<" fsnr = "<<fsnr<<" distsac = "<<distsac<<" famp = "<<famp<<" loaded = "<<loaded<<std::endl;
   }

   /* Initialize */
   void Initialize() {
      dist = azi1 = azi2 = snr = amp = grT = phT = -12345.;
   }

   /* compute/extract all available info */
   bool ExtractAll( float per ) {
      if( ! loaded ) return false;
		// compute distance and azimuths
		if( dist==-12345. || azi1==-12345. || azi2==-12345. )
			ComputeDisAzi();
      bool suc = Get_Amp_SNR( per );
      suc &= Get_GrT_PhT( per );
      if( hasAmpfile ) suc &= Get_Amp( per );
      return suc;
   }

   /* compute distance and azimuths */
	void ComputeDisAzi() {
		Path<float> path(src.lon, src.lat, sta.lon, sta.lat);
		dist = path.Dist(); azi1 = path.Azi1(); azi2=path.Azi2();
	}
   
	/* compute amplitude at period per from file famp */
	bool Get_Amp( float per ) {
		std::ifstream fin(famp);
		if( ! fin ) return false;
		float perlst[300], amplst[300], ftmp;
		int n = 0;
		for(std::string line; std::getline(fin, line); ) {
			if(sscanf(line.c_str(), "%f %f", &perlst[n], &amplst[n])!=2) continue;
			//amplst[n] *= 25;
			n++;
		}
		fin.close();
		InserSort(perlst,amplst,n);
		int i;
		for(i=0;i<n && perlst[i]<per; i++) {}
		if(i==0 || i==n) return false;
		amp = amplst[i-1]+(per-perlst[i-1])/(perlst[i]-perlst[i-1])*(amplst[i]-amplst[i-1]);
		return true;
	}

	/* compute amplitude and snr at period per from file fsnr */
	bool Get_Amp_SNR( float per ) {
		std::ifstream fin(fsnr.c_str());
		if( ! fin ) return false;
		float perlst[300], amplst[300], snrlst[300];
		int n = 0;
		for(std::string line; std::getline(fin, line); ) {
			if(sscanf(line.c_str(), "%f %f %f", &perlst[n], &amplst[n], &snrlst[n])!=3) continue;
			n++;
		}
		fin.close();
		InserSort2(perlst,amplst,snrlst,n);
		int i;
		for(i=0;i<n && perlst[i]<per; i++) {}
		if(i==0 || i==n) return false;
		amp = amplst[i-1]+(per-perlst[i-1])/(perlst[i]-perlst[i-1])*(amplst[i]-amplst[i-1]);
		snr = snrlst[i-1]+(per-perlst[i-1])/(perlst[i]-perlst[i-1])*(snrlst[i]-snrlst[i-1]);
		return true;
	}

	/* compute group and phase travel time at period per from file fdisp */
	bool Get_GrT_PhT( float per ) {
		std::ifstream fin(fdisp.c_str());
		if( ! fin ) return false;
		float perlst[300], grvlst[300], phvlst[300], ftmp;
		int itmp, n = 0;
		for(std::string line; std::getline(fin, line); ) {
			if(sscanf(line.c_str(), "%d %f %f %f %f", &itmp, &ftmp, &perlst[n], &grvlst[n], &phvlst[n])!=5) continue;
			n++;
		}
		fin.close();
		InserSort2(perlst,grvlst,phvlst,n);
		int i;
		for(i=0;i<n && perlst[i]<per; i++) {}
		if(i==0 || i==n) return false;
		grT = distsac / ( grvlst[i-1]+(per-perlst[i-1])/(perlst[i]-perlst[i-1])*(grvlst[i]-grvlst[i-1]) );
		phT = distsac / ( phvlst[i-1]+(per-perlst[i-1])/(perlst[i]-perlst[i-1])*(phvlst[i]-phvlst[i-1]) );
		return true;
	}

};



int main( int argc, char *argv[] ) {
	/* check input */
	if(argc!=5) {
		std::cerr<<"Usage: "<<argv[0]<<" [event.lst] [station.lst] [event-sta-fdisp-fsnr-distsac-famp list] [period]"<<std::endl;
		exit(-1);
	}

   /* read in event list */
   std::vector<PointC> evelst;
   std::ifstream fin(argv[1]);
   if( ! fin ) {
      std::cerr<<"Cannot read from file "<<argv[1]<<std::endl;
      exit(0);
   }
   for(std::string line; std::getline(fin, line); ) evelst.push_back( PointC(line.c_str()) );
   std::cerr<<evelst.size()<<" events read in."<<std::endl;
   if( evelst.size() == 0 ) exit(0);
   fin.close(); fin.clear();

   /* read in station list */
   std::vector<PointC> stalst;
   fin.open(argv[2]);
   if( ! fin ) {
      std::cerr<<"Cannot read from file "<<argv[2]<<std::endl;
      exit(0);
   }
   for(std::string line; std::getline(fin, line); ) stalst.push_back( PointC(line.c_str()) );
   std::cerr<<stalst.size()<<" stations read in."<<std::endl;
   if( stalst.size() == 0 ) exit(0);
   fin.close(); fin.clear();

   /* read in path list and fill in src/sta location from stalst and evelst */
   std::vector<PathRec> pathlst;
   fin.open(argv[3]);
   if( ! fin ) {
      std::cerr<<"Cannot read from file "<<argv[3]<<std::endl;
      exit(0);
   }
   for(std::string line; std::getline(fin, line); ) {
      PathRec pathtmp( stalst, evelst, line.c_str() );
      if( pathtmp.loaded ) pathlst.push_back( pathtmp );
   }
   std::cerr<<pathlst.size()<<" paths read in."<<std::endl;
   if( pathlst.size() == 0 ) exit(0);
   fin.close(); fin.clear();

   /* work on each path */
   float per = atof(argv[4]);
   char outname[30];
   sprintf(outname, "Disp_info_%.1fsec", per);
   std::ofstream fout(outname);
   fout<<pathlst.at(0).Format()<<std::endl;
   for(int i=0; i<pathlst.size(); i++) {
      pathlst.at(i).ExtractAll(per);
      fout<<pathlst.at(i)<<std::endl;
   }
   fout.close();
   std::cout<<outname<<std::endl;

   return 0;
}
