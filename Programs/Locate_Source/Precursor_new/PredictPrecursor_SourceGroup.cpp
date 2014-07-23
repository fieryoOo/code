/* this code takes 
 * (1) a SAC file list (just to be convenient. The 2 station locations stored in the SAC header will be used to predict the arrival time)
 * (2) a source group list (source_lon source_lat), and
 * (3) a group velocity map of the region of interests, and 
 * and predict the arrival time of the precursoring signal for each path */
#include "SacRec.h"
#include "Map.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <new>
#include <vector>
#include <stdexcept>

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
			char sacname[300];
			if( sscanf(line.c_str(), "%s", sacname)!= 1 )
				throw std::runtime_error( std::string("Error: Format error in sac list: ") + line );
			_saclst.push_back(sacname);
		}
		fin.close();
		_icur = 0;
	}

	bool PrepareNext() {
		if( _icur < 0 || 
			 _saclst.size() <= 0 ||
			 _icur >= _saclst.size() ) return false;

		SacRec sac( _saclst.at(_icur).c_str() );
		sac.LoadHD();
		_dist = sac.shd.dist;
		_azi = sac.shd.az;
		_baz = sac.shd.baz;
		float evlo = sac.shd.evlo<0. ? sac.shd.evlo+360. : sac.shd.evlo;
		float stlo = sac.shd.stlo<0. ? sac.shd.stlo+360. : sac.shd.stlo;
		_Psta1 = Point<float>( evlo, sac.shd.evla );
		_Psta2 = Point<float>( stlo, sac.shd.stla );
		_name = _saclst[_icur];
		_icur++;
		return true;
	}

	bool Rewind() {
		if( _icur < 0 || _saclst.size() <= 0 ) return false;
		_icur = 0;
		return true;
	}

	const char* Name() { return _name.c_str(); }
	float Dist() { return _dist; }
	float Azi() { return _azi; }
	float Baz() { return _baz; }
	Point<float> P1() { return _Psta1; }
	Point<float> P2() { return _Psta2; }

	size_t size() { return _saclst.size(); }

private:
	std::vector<std::string> _saclst;
	Point<float> _Psta1, _Psta2;
	std::string _name;
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

private:
	std::vector< Point<float> > _srclst;
	int _icur;
	Point<float> _P;
	static constexpr float NaN = -12345.;

};


#include <limits>
struct PathResult {
	PathResult() 
		: tsw(NaN)
		, tpmin( std::numeric_limits<float>::max() )
		, tpmax( std::numeric_limits<float>::lowest() ) {}

   friend std::ostream& operator << (std::ostream& o, PathResult a) { 
      o << a.tsw << "   " << a.tpmin << " " <<a.tpmax; 
      return o; 
   }

	float tsw, tpmin, tpmax;
	static constexpr float NaN = -12345.;
};

#define DEBUG
int main( int argc, char* argv[] ) 
{
	try {
		/* check for and read in input parameters */
		if( argc != 5 ) {
			std::cerr<<"Usage: "<<argv[0]<<" [SAC list] [source list (slon slat)] [Input model] [Out file]"<<std::endl;
			exit(-1);
		}
		char* mapname = argv[2];

		/* load in sac list */
		SACHandler saclst( argv[1] );

		/* load in source locations */
		SourceHandler srclst( argv[2] );


		/* loop through all sources and all sac files and compute t_precursors */
		std::vector<PathResult> results( saclst.size() );
		float perc1, perc2, lamda = 15.;
		while( srclst.PrepareNext() ) {
			/* define the map object */
			Map velmap( argv[3], srclst.P() );
			saclst.Rewind();
			int isac = 0;
			while( saclst.PrepareNext() ) {
				DataPoint<float> res1 = velmap.PathAverage_Reci( saclst.P1(), lamda, perc1 );
				DataPoint<float> res2 = velmap.PathAverage_Reci( saclst.P2(), lamda, perc2 );
				if( perc1 > 0.9 && perc2 > 0.9 ) {
					float t1 = res1.Dis() / res1.Data();
					float t2 = res2.Dis() / res2.Data();
					float tp = t2 - t1;
					PathResult& presult = results.at(isac);
					if( presult.tpmin > tp ) presult.tpmin = tp;
					if( presult.tpmax < tp ) presult.tpmax = tp;
				}
				isac++;
			}
		}

		/* loop through all sac files, compute t_SW and output */
		char* outname = argv[4];
		std::ofstream fout(outname);
		saclst.Rewind();
		int isac = 0;
		fout<<"sacf_name(1) intersta_dist(2) azimuth(3) back_azi(4)   T_surfacewave(5)   T_precursor_min(6) T_precursor_max(7)"<<std::endl;
		while( saclst.PrepareNext() ) {
			Map velmap( argv[3], saclst.P1() );
			DataPoint<float> res = velmap.PathAverage_Reci( saclst.P2(), lamda, perc1 );
			if( perc1 > 0.9 ) {
				PathResult& presult = results.at(isac);
				presult.tsw = res.Dis() / res.Data();
				fout<<saclst.Name()<<" "<<saclst.Dist()<<" "<<saclst.Azi()<<" "<<saclst.Baz()<<"   "<<presult<<std::endl;
			}
			isac++;
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
