/* this program takes 
 * (1) a list of SAC files (in the header of which 'EVLO EVLA STLO STLA' must present)
 * (2) a source group list (source_lon source_lat), and
 * (3) a group velocity map of the region of interests,
 * predict the arrival time window of the precursoring signal,
 * compute and output the rmsPrecursor/rmsNoise for each SAC file */
#include "SacRec.h"
#include "Map.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <new>
#include <vector>
#include <stdexcept>
#include <limits>

class SACHandler {

public:
	SACHandler( const char* fname ) 
		: _icur(NaN), _dist(NaN), _azi(NaN), _baz(NaN) {
		Load( fname );
	}

	void Load( const char* fname ) {
		_saclst.clear();
		std::ifstream fin(fname);
		if( !fin ) 
			throw std::runtime_error( std::string("Error: Cannot access sac list file ") + fname );
		for( std::string line; std::getline(fin, line); ) {
			char sacname[300], otherinfo[300];
			if( sscanf(line.c_str(), "%s %300[^\n]", sacname, otherinfo) < 1 )
				throw std::runtime_error( std::string("Error: Format error in sac list: ") + line );
			_saclst.push_back(sacname);
			_otherinfo.push_back(otherinfo);
		}
		fin.close();
		_icur = 0;
	}

	bool PrepareNext() {
		if( _icur < 0 || 
			 _saclst.size() <= 0 ||
			 _icur >= _saclst.size() ) return false;

		SacRec sac( _saclst.at(_icur) );
		sac.LoadHD();
		_dist = sac.shd.dist;
		_azi = sac.shd.az;
		_baz = sac.shd.baz;
		float evlo = sac.shd.evlo<0. ? sac.shd.evlo+360. : sac.shd.evlo;
		float stlo = sac.shd.stlo<0. ? sac.shd.stlo+360. : sac.shd.stlo;
		_Psta1 = Point<float>( evlo, sac.shd.evla );
		_Psta2 = Point<float>( stlo, sac.shd.stla );
		_name = _saclst[_icur];
		_oinfo = _otherinfo[_icur];
		_icur++;
		return true;
	}

	bool Rewind() {
		if( _icur < 0 || _saclst.size() <= 0 ) return false;
		_icur = 0;
		return true;
	}

	float ComputeCurrentSNR( float tpmin, float tpmax, float tsw ) {
		SacRec sac( _name );
		sac.Load();
		float rmsPrecursor, rmsNoise = std::numeric_limits<float>::max();
		sac.RMSAvg( tpmin, tpmax, rmsPrecursor );
		// find the quietest window to be the noise window
		int winstep = 500, winlen=800;
		int sign = -1;
		float tnmin = -tsw-300.-winlen, tnmax;
		while( true ) {
			try {
				float rmscur;
				sac.RMSAvg( tnmin, tnmax, rmscur );
				if( rmsNoise > rmscur ) rmsNoise = rmscur;
				tnmin += sign * winstep;
				tnmax = tnmin + winlen;
			} catch(std::exception& e) {
				if( tnmin < 0. ) {
					tnmin = tsw+300.;
					sign = 1;
				} else break;
			}
		}

		if( rmsNoise == std::numeric_limits<float>::max() ) {
			return -12345.;
		} else {
			float daynum = sac.shd.user0;
			if( daynum <= 0 || daynum != (int)daynum )
				throw std::runtime_error(" #day info not available in sac header (user0)!");
			//return 10.*rmsPrecursor / (rmsNoise*sqrt(daynum));
			return rmsPrecursor;
		}
	}

	const char* Name() { return _name.c_str(); }
	const char* OInfo() { return _oinfo.c_str(); }
	float Dist() { return _dist; }
	float Azi() { return _azi; }
	float Baz() { return _baz; }
	Point<float> P1() { return _Psta1; }
	Point<float> P2() { return _Psta2; }

	size_t size() { return _saclst.size(); }

private:
	std::vector<std::string> _saclst;
	std::vector<std::string> _otherinfo;
	Point<float> _Psta1, _Psta2;
	std::string _name, _oinfo;
	int _icur;
	float _dist, _azi, _baz;
	static constexpr float NaN = -12345.;

};


class SourceHandler {

public:
	SourceHandler( const char* fname ) 
		: _P(NaN), _icur(NaN) {
		Load( fname );
	}

	void Load( const char* fname ) {
		_srclst.clear();
		std::ifstream fin(fname);
		if( !fin ) 
			throw std::runtime_error( std::string("Error: Cannot access sac list file ") + fname );
		for( std::string line; std::getline(fin, line); ) {
			float slon, slat;
			if( sscanf(line.c_str(), "%f %f", &slon, &slat) != 2 )
				throw std::runtime_error( std::string("Error: Format error in source list: ") + line );
			if( slon < 0.) slon += 360.;
			_srclst.push_back( Point<float>(slon, slat) );
			//std::cerr<<*(_srclst.end()-1)<<std::endl;
		}
		fin.close();
		_icur = 0;
	}

	bool PrepareNext() {
		if( _icur < 0 || 
			 _srclst.size() <= 0 ||
			 _icur >= _srclst.size() ) return false;
		_P = Point<float>( _srclst[_icur].Lon(), _srclst[_icur].Lat() );
		_icur++;
		return true;
	}

	bool Rewind() {
		if( _icur < 0 || _srclst.size() <= 0 ) return false;
		_icur = 0;
		return true;
	}

	Point<float> P() { return _P;	}

	size_t size() { return _srclst.size(); }

private:
	std::vector< Point<float> > _srclst;
	int _icur;
	Point<float> _P;
	static constexpr float NaN = -12345.;

};


struct PathResult {

	bool validTPs() {
		//return ( tpmin<std::numeric_limits<float>::max() && tpmax>std::numeric_limits<float>::lowest() );
		return tpmin<tpmax;
	}

   friend std::ostream& operator << (std::ostream& o, PathResult a) { 
      o << a.tsw << "   " << a.tpmin << " " << a.tpmax << "   " << a.srDis1/a.count << " " << a.srDis2/a.count << " " << a.snr; 
      return o; 
   }

	int count;
	float srDis1 = 0, srDis2 = 0;
	float tsw = NaN, snr = NaN;
	float tpmin = std::numeric_limits<float>::max();
	float tpmax = std::numeric_limits<float>::lowest();
	static constexpr float NaN = -12345.;
};

#define DEBUG
int main( int argc, char* argv[] ) 
{
	try {
		/* check for and read in input parameters */
		if( argc != 5 && argc != 6 ) {
			std::cerr<<"Usage: "<<argv[0]<<" [SAC list] [source list (slon slat)] [Input model] [Out file] "
						<<"[(optional) Tpred outfile name (for all source points)]"<<std::endl;
			exit(-1);
		}
		char* mapname = argv[2];

		/* load in sac list */
		SACHandler saclst( argv[1] );
		std::cout<<saclst.size()<<" sac file(s) loaded."<<std::endl;

		/* load in source locations */
		SourceHandler srclst( argv[2] );
		std::cout<<srclst.size()<<" source(s) loaded."<<std::endl;


		/* loop through all sources and all sac files and compute t_precursors */
		bool writeout2 = false;
		std::ofstream fout2;
		if( argc == 6 ) {
			fout2.open(argv[5]);
			writeout2 = true;
		}
		std::vector<PathResult> results( saclst.size() );
		float perc1, perc2, lamda = 1.;
		for( int isrc=0; srclst.PrepareNext(); isrc++ ) {
			std::cout<<"   Working on source #"<<isrc<<std::endl;
			/* define the map object */
			Map velmap( argv[3], srclst.P() );
			saclst.Rewind();
			for( int isac=0; saclst.PrepareNext(); isac++ ) {
				DataPoint<float> res1 = velmap.PathAverage_Reci( saclst.P1(), lamda, perc1 );
				DataPoint<float> res2 = velmap.PathAverage_Reci( saclst.P2(), lamda, perc2 );
				if( perc1 > 0.9 && perc2 > 0.9 ) {
					float t1 = res1.Dis() / res1.Data();
					float t2 = res2.Dis() / res2.Data();
					float tp = t2 - t1;
					PathResult& presult = results.at(isac);
					if( presult.tpmin > tp ) presult.tpmin = tp;
					if( presult.tpmax < tp ) presult.tpmax = tp;
					presult.srDis1 += res1.Dis();
					presult.srDis2 += res2.Dis();
					presult.count++;
					if( writeout2) fout2<<srclst.P()<<" "<<tp<<"   "<<saclst.Name()<<std::endl;
				}
			}
		}
		fout2.close();

		/* loop through all sac files, compute t_SW, snrPrecursor, and output */
		char* outname = argv[4];
		std::ofstream fout(outname);
		saclst.Rewind();
		fout<<"sacf_name(1) inters_dis(2) azimuth(3) back_azi(4)   T_sw(5)   T_p_min(6) T_p_max(7)   s-r_dis1(8) s-r_dis2(9) SNRPrecursor(10) OtherInfo(11 on)"<<std::endl;
		//int Tsafe = 60;
		for( int isac=0; saclst.PrepareNext(); isac++ ) {
			Map velmap( argv[3], saclst.P1() );
			DataPoint<float> res = velmap.PathAverage_Reci( saclst.P2(), lamda, perc1 );
			if( perc1 > 0.9 ) {
				PathResult& presult = results.at(isac);
				if( ! presult.validTPs() ) continue;
				presult.tsw = res.Dis() / res.Data();
				//float TPmax = presult.tsw - Tsafe;
				//if( presult.tpmin<-TPmax || presult.tpmax>TPmax ) continue;
				presult.snr = saclst.ComputeCurrentSNR( presult.tpmin, presult.tpmax, presult.tsw );
				if( presult.snr > 0. )
					fout<<saclst.Name()<<" "<<res.Dis()<<" "<<saclst.Azi()<<" "<<saclst.Baz()<<"   "<<presult<<"   "<<saclst.OInfo()<<std::endl;
				else
					std::cerr<<"Warning(main): failed when computing SNR for sacfile "<<saclst.Name()<<std::endl;
			}
		}
		fout.close();

	} catch(std::exception& e) {
		std::cerr<<e.what()<<std::endl;
		return -1;
	} catch(...) {
		std::cerr<<"Unknown Exception!"<<std::endl;
		return -2;
	}

	return 0;
}
