#include "Map.h"
#include "DisAzi.h"
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

template <typename T>
class ppair {
public:
	int pos = -1;	// complementary position info
	//ppair() {}
	ppair( const typename std::vector<T>::iterator p1, const typename std::vector<T>::iterator p2 )
		: p1(p1), p2(p2) {}
	ppair( const typename std::vector<T>::iterator p1, const size_t size )
		: p1(p1), p2(p1+size) {}
	size_t size() const { return p2-p1; }
	inline const typename std::vector<T>::iterator& begin() const { return p1; }
	inline		 typename std::vector<T>::iterator& begin() { return p1; }
	inline const typename std::vector<T>::iterator& end() const { return p2; }
	inline		 typename std::vector<T>::iterator& end() { return p2; }
private:
	typename std::vector<T>::iterator p1, p2;
};

struct Map::Mimpl {

	Mimpl() {}

	Mimpl( const Mimpl& m2 )
		: lonmin(m2.lonmin), lonmax(m2.lonmax), latmin(m2.latmin), latmax(m2.latmax)
		, dis_lat1D(m2.dis_lat1D), dis_lon1D(m2.dis_lon1D), dataV(m2.dataV), dataavg(m2.dataavg)
		, grd1_lon(m2.grd1_lon), grd1_lat(m2.grd1_lat), dataM1(m2.dataM1)
		, grd2_lon(m2.grd2_lon), grd2_lat(m2.grd2_lat), dataM2(m2.dataM2), isReg(m2.isReg) {
		ApplyHash();
	}

	Mimpl( Mimpl&& m2 )
		: lonmin(m2.lonmin), lonmax(m2.lonmax), latmin(m2.latmin), latmax(m2.latmax)
		, dis_lat1D(m2.dis_lat1D), dis_lon1D(std::move(m2.dis_lon1D)), dataV(std::move(m2.dataV)), dataavg(std::move(m2.dataavg))
		, grd1_lon(m2.grd1_lon), grd1_lat(m2.grd1_lat), dataM1(std::move(m2.dataM1))
		, grd2_lon(m2.grd2_lon), grd2_lat(m2.grd2_lat), dataM2(std::move(m2.dataM2)), isReg(m2.isReg) {
		// this is likely not necessary when moving unless the compiler writer went nuts
		ApplyHash();	// but just to be safe...
	}

	Mimpl& operator= ( const Mimpl& m2 ) {
		lonmin=m2.lonmin; lonmax=m2.lonmax; latmin=m2.latmin; latmax=m2.latmax;
		dis_lat1D=m2.dis_lat1D; dis_lon1D=m2.dis_lon1D; dataV=m2.dataV; dataavg=m2.dataavg;
		grd1_lon=m2.grd1_lon; grd1_lat=m2.grd1_lat; dataM1=m2.dataM1;
		grd2_lon=m2.grd2_lon; grd2_lat=m2.grd2_lat; dataM2=m2.dataM2; isReg=m2.isReg;
		ApplyHash();
		return *this;
	}

	Mimpl& operator= ( Mimpl&& m2 ) {
		lonmin=m2.lonmin; lonmax=m2.lonmax; latmin=m2.latmin; latmax=m2.latmax;
		dis_lat1D=m2.dis_lat1D; dis_lon1D=m2.dis_lon1D; dataV=std::move(m2.dataV); dataavg=std::move(m2.dataavg);
		grd1_lon=m2.grd1_lon; grd1_lat=m2.grd1_lat; dataM1=std::move(m2.dataM1);
		grd2_lon=m2.grd2_lon; grd2_lat=m2.grd2_lat; dataM2=std::move(m2.dataM2); isReg=m2.isReg;
		// this is likely not necessary when moving unless the compiler writer went nuts
		ApplyHash();	// but just to be safe...
		return *this;
	}


	float lonmin = NaN, lonmax = NaN;
	float latmin = NaN, latmax = NaN;

	float dis_lat1D;
	std::vector<float> dis_lon1D;

	// data vector
	std::vector< DataPoint<float> > dataV;
	std::vector< DataPoint<float> > dataavg;

	// matrix (of iterators) 1
	//Array2D< std::vector< DataPoint<float> > > dataM;
	float grd1_lon = 1., grd1_lat = 1.;
	Array2D< ppair< DataPoint<float> > > dataM1;
	// matrix (of iterators) 2
	bool isReg = false;
	float grd2_lon, grd2_lat;
	Array2D< ppair< DataPoint<float> > > dataM2;


	inline int ilon( float xloc ) const { float xgrd = isReg ? grd2_lon : grd1_lon; return (int)floor((xloc-lonmin) / xgrd + 0.5); }
	inline int ilat( float yloc ) const { float ygrd = isReg ? grd2_lat : grd1_lat; return (int)floor((yloc-latmin) / ygrd + 0.5); }
	inline int ilon_floor( float xloc ) const { float xgrd = isReg ? grd2_lon : grd1_lon; return (int)floor(roundoff(xloc-lonmin) / xgrd); }
	inline int ilat_floor( float yloc ) const { float ygrd = isReg ? grd2_lat : grd1_lat; return (int)floor(roundoff(yloc-latmin) / ygrd); }
	inline int ilon_ceil( float xloc ) const { float xgrd = isReg ? grd2_lon : grd1_lon; return (int)ceil(roundoff(xloc-lonmin) / xgrd); }
	inline int ilat_ceil( float yloc ) const { float ygrd = isReg ? grd2_lat : grd1_lat; return (int)ceil(roundoff(yloc-latmin) / ygrd); }


	void ReadData( const std::string& fname ) {
		std::ifstream fin(fname.c_str());
		if( ! fin )
			throw ErrorM::BadFile(FuncName, "read from "+fname);
		float lon, lat, data;
		//sscanf(line.c_str(), "%f %f %f", &lon, &lat, &data);
		//fin.seekg(0); // rewind
		/* load in all data */
		dataV.clear(); dataavg.resize(1, 0.);
		for(std::string line; std::getline(fin, line); ) {
			float lon, lat, data;
			std::stringstream ss(line);
			if( ! ( ss >> lon >> lat >> data) ) {
				//std::cerr<<"Warning(Map::Load): wrong format detected in file "<<fname<<" ("<<line<<") "<<std::endl;
				continue;
			}
			if( data != data ) continue;
			if( lon < 0. ) lon += 360.;
			//float dis = Path<float>(srcin, Point<float>(lon,lat)).Dist();
			dataV.push_back( DataPoint<float>(lon, lat, data) );
			dataavg[0].data += data;
		}
		fin.close();
		dataavg[0].data /= dataV.size();
		if( dataV.size() == 0 )
			throw ErrorM::BadFile(FuncName, "empty file "+fname);

		//std::cout<<"### "<<dataV.size()<<" map points loaded. ###"<<std::endl;
		/* /map boundaries */
		CompMapBoundaries();
	}

	void Hash() {
		if( dataV.size() == 0 ) return;
		HashM1();
		isReg = HashM2();
		ApplyHash();
	}

	/* matrix #1: store by blocks */
	void HashM1() {
		/* grid/hash size */
		//grd1_lon = grd1_lat = sqrt( (lonmax-lonmin) * (latmax-latmin) / sqrt((float)dataV.size()) );
		//grd1_lon = grdlon; grd1_lat = grdlat;
		int nlon = (int)ceil( (lonmax-lonmin) / grd1_lon ) + 1;
		int nlat = (int)ceil( (latmax-latmin) / grd1_lat ) + 1;
		Array2D< std::vector< DataPoint<float> > > dataM_hash;
		dataM_hash.clear(); dataM_hash.resize(nlon, nlat);
		int size = 0; dataavg.resize(1, 0.);
		Point<float> BL(lonmin, latmin);
		Point<float> TR(lonmax, latmax);
		for( const auto& dpcur : dataV ) {
			if( ! dpcur.isWithin(BL, TR) ) continue;
			int ilon = (int)floor((dpcur.Lon()-lonmin) / grd1_lon + 0.5);
			int ilat = (int)floor((dpcur.Lat()-latmin) / grd1_lat + 0.5);
			dataM_hash(ilon, ilat).push_back(dpcur);
			dataavg[0].data += dpcur.data;
			size++;
		}
		dataavg[0].data /= size;
		/* now rearrange all DataPoints into dataV and store iterators in dataM */
		dataV.clear(); dataV.reserve(size*2+1000);
		dataM1.resize(nlon, nlat, ppair< DataPoint<float> >( dataavg.begin(), 0 ));
		for( size_t i=0; i<dataM_hash.NumRows(); i++ ) {
			for( size_t j=0; j<dataM_hash.NumCols(); j++ ) {
				auto& dpV = dataM_hash(i, j);
				/* undefined behavier?
				auto p1 = dataV.end();
				dataV.insert( dataV.end(), std::make_move_iterator(dpV.begin()), 
													std::make_move_iterator(dpV.end()) );
				auto p2 = dataV.end();
				*/
				//size_t old_size = dataV.size(), inser_size = dpV.size();
				dataM1(i, j).pos = dataV.size();
				std::copy( std::make_move_iterator(dpV.begin()), 
							  std::make_move_iterator(dpV.end()), std::back_inserter(dataV) );
				//dataM1(i, j) = ppair< DataPoint<float> >(dataV.begin() + old_size, inser_size);
			}
		}
		/* compute distance of 1 degree in lon/lat */
		//dis_lat1D = 111.;
		float clon = 0.5 * (lonmin + lonmax);
		dis_lat1D = Path<float>(clon, 0., clon, 1.).Dist();
		//dis_lon1D.resize( dataM.NumCols() );
		dis_lon1D.clear();
		float latcur = latmin;
		while( latcur<=90. ) {
			dis_lon1D.push_back( Path<float>(0., latcur, 1., latcur).Dist() );
			latcur += grd1_lat;
		}
	}

	/* matrix #2: store by singles */
	bool HashM2() {
		if( ! CompMapGrids() ) return false;
		int nlon = (int)ceil( (lonmax-lonmin) / grd2_lon ) + 1;
		int nlat = (int)ceil( (latmax-latmin) / grd2_lat ) + 1;
		dataM2.clear(); 
		dataM2.resize(nlon, nlat, ppair< DataPoint<float> >( dataavg.begin(), 1 ) );
		//for( auto iter=dataV.begin(); iter<dataV.end(); iter++ ) {
		for(int i=0; i<dataV.size(); i++) {
			const auto& dp = dataV[i];
			float ilon = (int)floor( (dp.lon-lonmin) / grd2_lon + 0.5 );
			float ilat = (int)floor( (dp.lat-latmin) / grd2_lat + 0.5 );
         //dataM2(ilon, ilat) = ppair< DataPoint<float> >(dataV.begin()+i, 1);
			dataM2(ilon, ilat).pos = i;
      }
		if( (float)dataV.size() / dataM2.Size() < 0.9 ) {
			dataM2.clear();
			return false;
		}
		return true;
	}

	void ApplyHash() {
		// apply hash on dataM1
		int i;
		for( i=0; i<dataM1.Size()-1; i++ ) {
			auto& pp = dataM1[i];
			if( pp.pos == -1 ) continue;
			pp = ppair< DataPoint<float> >( dataV.begin()+pp.pos, dataV.begin()+dataM1[i+1].pos );
		}
		dataM1[i] = ppair< DataPoint<float> >( dataV.begin()+dataM1[i].pos, dataV.end() );

		// apply hash on dataM2
		if( ! isReg ) return;
		for( auto& pp : dataM2 ) {
			if( pp.pos == -1. ) continue;
			pp = ppair< DataPoint<float> >( dataV.begin()+pp.pos, 1 );
		}
	}

   // compute map grids and boundaries
	inline float roundoff( float val ) const { return (int)floor(val*10000+0.5)/10000.; }

	void CompMapBoundaries( const float lomin = -99999., const float lomax = 99999.,
									const float lamin = -99999., const float lamax = 99999. ) {
		Point<float> BL(lomin, lamin);
		Point<float> TR(lomax, lamax);
		bool notinit = true;
		for(size_t i=0; i<dataV.size(); i++) {
			const auto& dpcur = dataV[i];
			if( ! dpcur.isWithin(BL, TR) ) continue;
			if( notinit ) {
				lonmin = lonmax = dpcur.Lon();
				latmin = latmax = dpcur.Lat();
				notinit = false;
			} else {
				if( lonmax < dpcur.Lon() ) {
					lonmax = dpcur.Lon();
				} else if( lonmin > dpcur.Lon() ) {
					lonmin = dpcur.Lon();
				}
				if( latmax < dpcur.Lat() ) {
					latmax = dpcur.Lat();
				} else if( latmin > dpcur.Lat() ) {
					latmin = dpcur.Lat();
				}
			}
		}
	}

   bool CompMapGrids() {
      // compute x grid
		// sort by lon
      auto dataVtmp = dataV;
      std::sort( dataVtmp.begin(), dataVtmp.end(), [](const DataPoint<float>& n1, const DataPoint<float>& n2) {
            return n1.lon < n2.lon;
      } );
		// remove redundances
      auto last = std::unique(dataVtmp.begin(), dataVtmp.end(), [&](const DataPoint<float>& n1, const DataPoint<float>& n2) {
         return roundoff(n1.lon-n2.lon) == 0;
      } );
      dataVtmp.erase(last, dataVtmp.end());
		// search for smallest diff
      //xmin = dataVtmp.front().x;
      //xmax = dataVtmp.back().x;
      grd2_lon = roundoff(dataVtmp[1].lon - dataVtmp[0].lon);
      for(int i=2; i<dataVtmp.size(); i++) {
         float xgrdc = roundoff(dataVtmp[i].lon-dataVtmp[i-1].lon);
         if( grd2_lon > xgrdc ) grd2_lon = xgrdc;
      }
		// check for irregular grids
      for(int i=2; i<dataVtmp.size(); i++) {
         float xgrdc = roundoff(dataVtmp[i].lon-dataVtmp[i-1].lon);
			float mulf = xgrdc / grd2_lon;
			int muli = (int)floor(mulf + 0.5);
			if( roundoff(mulf-muli)!=0 ) return false;
      }
      // compute y grid
		// sort by lon
      dataVtmp = dataV;
      std::sort( dataVtmp.begin(), dataVtmp.end(), [](const DataPoint<float>& n1, const DataPoint<float>& n2) {
            return n1.lat < n2.lat;
      } );
		// remove redundances
      last = std::unique(dataVtmp.begin(), dataVtmp.end(), [&](const DataPoint<float>& n1, const DataPoint<float>& n2) {
         return roundoff(n1.lat-n2.lat) == 0;
      } );
      dataVtmp.erase(last, dataVtmp.end());
		// search for smallest diff
      //ymin = dataVtmp.front().y;
      //ymax = dataVtmp.back().y;
      grd2_lat = roundoff(dataVtmp[1].lat - dataVtmp[0].lat);
      for(int i=2; i<dataVtmp.size(); i++) {
         float ygrdc = roundoff(dataVtmp[i].lat-dataVtmp[i-1].lat);
         if( grd2_lat > ygrdc ) grd2_lat = ygrdc;
      }
		// check for irregular grids
      for(int i=2; i<dataVtmp.size(); i++) {
         float ygrdc = roundoff(dataVtmp[i].lat-dataVtmp[i-1].lat);
			float mulf = ygrdc / grd2_lat;
			int muli = (int)floor(mulf + 0.5);
			if( roundoff(mulf-muli)!=0 ) return false;
      }

		return true;
   }

	inline float estimate_dist(const Point<float>& p1, const Point<float>& p2) {
		float lon1 = p1.Lon(), lat1 = p1.Lat();
		if( lon1 < 0. ) lon1 += 360.;
		float lon2 = p2.Lon(), lat2 = p2.Lat();
		if( lon2 < 0. ) lon2 += 360.;

		size_t idlmax = dis_lon1D.size() - 1;
		int ilatmid = (int)floor( ( (lat1+lat2) * 0.5 - latmin ) / grd1_lat + 0.5 );
		if( ilatmid > idlmax ) ilatmid = idlmax;
		else if( ilatmid < 0 ) ilatmid = 0;

		float dis_lon = fabs(lon1-lon2);
		if( dis_lon > 180. ) dis_lon = 360. - dis_lon;
		dis_lon *= dis_lon1D[ilatmid];
		float dis_lat = (lat1-lat2) * dis_lat1D;
/*
std::cerr<<fabs(lon1-lon2)<<" "<<(lat1-lat2)<<"\n";
std::cerr<<dis_lon1D[ilatmid]<<" "<<dis_lat1D<<"\n";
std::cerr<<p1<<" "<<p2<<"   "<<dis_lon<<" "<<dis_lat<<"\n";
*/
		return sqrt(dis_lon*dis_lon + dis_lat*dis_lat);
	}

	float Interp4( const DataPoint<float>& dp1, const DataPoint<float>& dp2,
						const DataPoint<float>& dp3, const DataPoint<float>& dp4, const Point<float>& P ) {
		float w1 = 1. / (0.01+estimate_dist( dp1, P )), dat1 = dp1.data;
		float w2 = 1. / (0.01+estimate_dist( dp2, P )), dat2 = dp2.data;
		float w3 = 1. / (0.01+estimate_dist( dp3, P )), dat3 = dp3.data;
		float w4 = 1. / (0.01+estimate_dist( dp4, P )), dat4 = dp4.data;
		return (w1*dat1 + w2*dat2 + w3*dat3 + w4*dat4 ) / (w1+w2+w3+w4);
	}

	float Interpolate( const Point<float>& P ) {
		const auto& dataM = isReg ? dataM2 : dataM1;
		int ilon_l = ilon_floor( P.lon ), ilon_u = ilon_ceil( P.lon );
		int ilat_l = ilat_floor( P.lat ), ilat_u = ilat_ceil( P.lat );
		if( ilon_l<0 || ilon_u>=dataM.NumRows() || ilat_l<0 || ilat_u>=dataM.NumCols() )
			return dataavg[0].data;
		if( isReg ) {
			return Interp4( *(dataM(ilon_l, ilat_l).begin()), *(dataM(ilon_l, ilat_u).begin()),
								 *(dataM(ilon_u, ilat_l).begin()), *(dataM(ilon_u, ilat_u).begin()), P );
		} else {
			float dismin = 99999.;
			float valbest = dataavg[0].data;
			for( int ilon=ilon_l; ilon<=ilon_u; ilon++ )
				for( int ilat=ilat_l; ilat<=ilat_u; ilat++ )
					for( const auto& dp : dataM(ilon,ilat) ) {
						float dis = estimate_dist(dp, P);
						if( dismin > dis ) {
							dismin = dis; valbest = dp.data;
						}
					}
			return valbest;
			/* too slow!
			std::vector< std::vector< DataPoint<float> >::const_iterator > iV;
			for( int ilon=ilon_l; ilon<=ilon_u; ilon++ )
				for( int ilat=ilat_l; ilat<=ilat_u; ilat++ )
					for( auto idp=dataM(ilon,ilat).begin(); idp<dataM(ilon,ilat).end(); idp++ ) {
						iV.push_back(idp);
					}
			if( iV.size() < 4 ) {
				//std::cerr<<"Warning(Map::Interpolate): no enough constraints at "<<P<<". Global average used"<<std::endl;
				return dataavg[0].data;
			}
			std::sort( iV.begin(), iV.end(), [&]( std::vector< DataPoint<float> >::const_iterator i1, 
															  std::vector< DataPoint<float> >::const_iterator i2 ) {
					return estimate_dist(*i1, P) < estimate_dist(*i2, P);
			} );
			return Interp4( *(iV[0]), *(iV[1]), *(iV[2]), *(iV[3]), P );
			*/
		}
	}

};


/* -------------- con/destructors and assignment operators ----------------- */
Map::Map( const float grdlon, const float grdlat ) 
	: pimplM( new Mimpl() ) {
	pimplM->grd1_lon = grdlon;
	pimplM->grd1_lat = grdlat;
}

Map::Map( const std::string& inname, const float grdlon, const float grdlat )
	: pimplM( new Mimpl() ), fname(inname) {
   pimplM->grd1_lon = grdlon;
   pimplM->grd1_lat = grdlat;
	Load( fname );
}
//	: Map( inname, Point<float>(), grdlon, grdlat) {}

Map::Map( const std::string& inname, const Point<float>& srcin, const float grdlon, const float grdlat ) 
	: pimplM( new Mimpl() ), fname(inname), src(srcin) {
	pimplM->grd1_lon = grdlon;
	pimplM->grd1_lat = grdlat;
	Load( fname );
	SetSource( src );
}

Map::Map( const Map& mp_other ) 
	: pimplM( new Mimpl(*(mp_other.pimplM)) )
	, fname(mp_other.fname), src(mp_other.src) {
	//pimplM->Hash(); // ApplyHash called inside Mimpl constructors
}

Map::Map( Map&& mp_other ) 
	: pimplM( std::move(mp_other.pimplM) ) 
	, fname(mp_other.fname), src(mp_other.src) { 
	//pimplM->Hash(); // ApplyHash called inside Mimpl constructors
}

Map& Map::operator= ( const Map& mp_other ) {
	pimplM.reset( new Mimpl(*(mp_other.pimplM)) );
	fname = mp_other.fname;
	src = mp_other.src;
	//pimplM->Hash(); // ApplyHash called inside Mimpl constructors
	return *this;
}

Map& Map::operator= ( Map&& mp_other ) {
	pimplM = std::move(mp_other.pimplM);
	fname = std::move(mp_other.fname);
	src = std::move(mp_other.src);
	//pimplM->Hash(); // ApplyHash called inside Mimpl constructors
	return *this;
}

Map::~Map() {}



/* ----- map boundaries ----- */
float Map::LonMin() const { return pimplM->lonmin; }
float Map::LonMax() const { return pimplM->lonmax; }
float Map::LatMin() const { return pimplM->latmin; }
float Map::LatMax() const { return pimplM->latmax; }

bool Map::isReg() const { return pimplM->isReg; }

/* ------------ IO and resets ------------ */
void Map::Load( const std::string& fnamein ) {
	fname = fnamein;
	/* check/autodefine src location
	if( src.Lon() == NaN || src.Lat() == NaN ) {
		std::ifstream fin( fnamein );
		std::string line;
		float lon, lat, data;
		if( fin && std::getline(fin, line) ) {
			sscanf(line.c_str(), "%f %f %f", &lon, &lat, &data);
		} else {
			throw ErrorM::BadFile(FuncName, "read from "+fnamein);
		}
		lon = lon>90. ? lon-90. : lon+270.;
		lat = lat>0 ? lat-90. : lat+90.;
		src = Point<float>(lon, lat);
	} */
	// open/check the file
	pimplM->ReadData( fname );
	pimplM->Hash();
}

/* ------------ set source location ------------ */
void Map::SetSource( const Point<float>& srcin ) {
	src = srcin;
	for( auto& dp : pimplM->dataV ) {
		dp.dis = Path<float>( src, Point<float>(dp.lon, dp.lat) ).Dist();
	}
}

/* --- clip the map around the source location (to speed up the average methods) --- */
void Map::Clip( const float lonmin, const float lonmax, const float latmin, const float latmax ) {
//std::cerr<<"Map::Clip 1:  "<<*(pimplM->dataM1(0,0).begin())<<" "<<*(pimplM->dataM2(0,0).begin())<<" "<<pimplM->dataM1.NumRows()<<" "<<pimplM->dataM1.NumCols()<<" "<<*this<<" "<<pimplM->dataV.size()<<" "<<pimplM->dataV.data()<<std::endl;
	// reset bounds to within the given region
	pimplM->CompMapBoundaries( lonmin, lonmax, latmin, latmax );
	// and re-hash
	pimplM->Hash();
}

size_t Map::size() const { return pimplM->dataV.size(); }

/* ------------ compute number of points near the given location ------------ */
float Map::NumberOfPoints(Point<float> rec, const float xhdis, const float yhdis, float& loneff, float& lateff) const {
	/* references */
	bool isReg = pimplM->isReg;
	const auto& dataM = isReg ? pimplM->dataM2 : pimplM->dataM1;
	float lonmin = pimplM->lonmin, latmin = pimplM->latmin;
	float grdlon = isReg ? pimplM->grd1_lon : pimplM->grd2_lon;
	float grdlat = isReg ? pimplM->grd1_lat : pimplM->grd2_lat;

   // requested boundaries
   Point<float> BL(rec.lon-xhdis, rec.lat-yhdis);
   Point<float> TR(rec.lon+xhdis, rec.lat+yhdis);

	int Npoints = dataM.Size();
	loneff = lateff = 0;
	if( isReg ) {
		// define computation area
      int rowmin = std::max( 0, pimplM->ilon_ceil(BL.lon) );
      int rowmax = std::min( dataM.NumRows(), pimplM->ilon_floor(TR.lon)+1 );
      int colmin = std::max( 0, pimplM->ilat_ceil(BL.lat) );
      int colmax = std::min( dataM.NumCols(), pimplM->ilat_floor(TR.lat)+1 );
		// compute total number of points in given area
		for(int irow=rowmin; irow<rowmax; irow++)
			for(int icol=colmin; icol<colmax; icol++) {
				const auto& dpcur = *(dataM(irow, icol).begin());
				loneff += dpcur.Lon(); lateff += dpcur.Lat();
			}
	} else {
		// define computation area
		int rowmin = std::max( 0, pimplM->ilon(BL.lon) );
		int rowmax = std::min( dataM.NumRows(), pimplM->ilon(TR.lon) + 1 );
		int colmin = std::max( 0, pimplM->ilat(BL.lat) );
		int colmax = std::min( dataM.NumCols(), pimplM->ilat(TR.lat) + 1 );
		// compute total number of points in given area
		Npoints = 0;
		for(int irow=rowmin; irow<rowmax; irow++)
			for(int icol=colmin; icol<colmax; icol++)
				for( const auto& dpcur : dataM(irow, icol) ) {
					if( ! dpcur.isWithin(BL, TR) ) continue;
					loneff += dpcur.Lon(); lateff += dpcur.Lat();
					Npoints++;
				}
	}

	if( Npoints > 0 ) {
		loneff /= Npoints;
		lateff /= Npoints;
	} else {
		loneff = lateff = NaN;
	}

	return Npoints;
}

/* ------------ compute average value on the point rec ------------ */
float Map::PointAverage(Point<float> rec, float hdis, float& weit) {
	/* references */
	const auto& dataM = pimplM->dataM1;
	float lonmin = pimplM->lonmin, latmin = pimplM->latmin;
	float grdlon = pimplM->grd1_lon, grdlat = pimplM->grd1_lat;;

	// correct lon
	rec.correctLon();

	// define computation area
	float dismax = hdis * 2.5, dismax_s = dismax*dismax;
	size_t idlmax = pimplM->dis_lon1D.size() - 1;
	int ilatrec = (int)floor((rec.Lat()-latmin) / grdlat + 0.5);
	if( ilatrec > idlmax ) ilatrec = idlmax;
	else if( ilatrec < 0 ) ilatrec = 0;
	float dis_lon1D = pimplM->dis_lon1D[ilatrec], dis_lat1D = pimplM->dis_lat1D;
	//calc_dist( rec.Lat(), rec.Lon(), rec.Lat(), rec.Lon()+1., &dis_lon1D );
	//calc_dist( rec.Lat(), rec.Lon(), rec.Lat()+1., rec.Lon(), &dis_lat1D );
	float Rlon = dismax / dis_lon1D, Rlat = dismax / dis_lat1D;
	int rowmin = (int)floor((rec.Lon()-Rlon - lonmin) / grdlon + 0.5);
	if( rowmin < 0 ) rowmin = 0;
	int rowmax = (int)floor((rec.Lon()+Rlon - lonmin) / grdlon + 0.5) + 1;
	if( rowmax > dataM.NumRows() ) rowmax = dataM.NumRows();
	int colmin = (int)floor((rec.Lat()-Rlat - latmin) / grdlat + 0.5);
	if( colmin < 0 ) colmin = 0;
	int colmax = (int)floor((rec.Lat()+Rlat - latmin) / grdlat + 0.5) + 1;
	if( colmax > dataM.NumCols() ) colmax = dataM.NumCols();

	weit = 0.;
	float alpha = -0.5 / (hdis*hdis), datasum = 0.;
	for(int irow=rowmin; irow<rowmax; irow++) {
		for(int icol=colmin; icol<colmax; icol++) {
			for( const auto& dpcur : dataM(irow, icol) ) {
				if( dpcur.Data() == NaN ) continue;
				//distance from dpcur to DPrec.
				float xdis = ( rec.Lon() - dpcur.Lon() ) * dis_lon1D;
				float ydis = ( rec.Lat() - dpcur.Lat() ) * dis_lat1D;
				float disc_s = xdis*xdis + ydis*ydis;
				//std::cerr<<dpcur.Data()<<" "<<disc_s<<"  "<<dpcur.Lon()<<" "<<dpcur.Lat()<<"  "<<rec.Lon()<<" "<<rec.Lat()<<std::endl;
				if( disc_s > dismax_s ) continue;
				float weight = exp( alpha * disc_s );
				weit += weight;
				datasum += (dpcur.Data() * weight);
			}
		}
	}
	if(weit==0.) datasum = NaN;
	else datasum /= weit;

	return datasum;

}

template<class Functor>
void Map::TraceGCP( Point<float> Psrc, const Point<float>& Prec, float dis_step, const Functor& func ) {
	// check source/receiver locations
	Point<float> BL( pimplM->lonmin, pimplM->latmin);
	Point<float> TR( pimplM->lonmax, pimplM->latmax);
	/*
	if( ! (Psrc.isWithin(BL, TR) && Prec.isWithin(BL, TR)) ) {
		std::cerr<<"Warning(Map::TraceGCP): src/rec "<<Psrc<<"/"<<Prec<<" out of bounds. Global average used"<<std::endl;
		func(pimplM->dataavg[0].data);
		return;
	}
	*/

	// compute tracing params
	Path<float> path( Psrc, Prec );
	float dis = path.Dist(), azi1 = path.Azi1(), azi2 = path.Azi2();
	int n = std::max(npts_min, (int)float( dis / dis_step + 0.5 ));
	dis_step = dis / n;

	// trace
	//Psrc.correctLon();
	//float val = pimplM->Interpolate(Psrc); func(val);
	//std::cerr<<Psrc<<" "<<val<<std::endl;
	func( pimplM->Interpolate(Psrc) );
	for(int i=0; i<n; i++) {
		path.reset( Psrc, dis_step, azi1 );
		Psrc = path.P2();
		azi1 = path.Azi2();
		func( pimplM->Interpolate(Psrc) );
		//float val = pimplM->Interpolate(Psrc); func(val);
		//std::cerr<<Psrc<<" "<<val<<std::endl;
	}
}

/* ------------ compute average value along the path src-rec ------------ */
DataPoint<float> Map::PathAverage(Point<float> rec, float& perc, const float lambda, const bool acc) {
	// check source
	if( src == Point<float>() )
		throw ErrorM::BadParam(FuncName, "invalid src location");

	float dis = Path<float>(src, rec).Dist();
	if( lambda < 0. ) {
		throw ErrorM::BadParam(FuncName, "negative lambda");
	} else if ( lambda == 0. ) { // call direct ray tracing if lambda == 0
		//if( !pimplM->isReg )
		//	throw ErrorM::BadParam(FuncName, "lambda cannot be 0 when map is irregular");
		int N = 0; float avg = 0.;
		TraceGCP(src, rec, trace_step, [&](float val) { avg += val; N++; } );
		avg /= N; perc = 1.;
		return DataPoint<float>(rec, avg, dis);
	}

	// references
	const auto& dataM = pimplM->dataM1;
	float lonmin = pimplM->lonmin, latmin = pimplM->latmin;
	float grd_lon = pimplM->grd1_lon, grd_lat = pimplM->grd1_lat;

	/* rec parameters */
	size_t idlmax = pimplM->dis_lon1D.size() - 1;
	int ilatmid = (int)floor( ( (rec.Lat()+src.Lat()) * 0.5 - latmin ) / grd_lat + 0.5 );
	if( ilatmid > idlmax ) ilatmid = idlmax;
	else if( ilatmid < 0 ) ilatmid = 0;
	float dis_lon1D = pimplM->dis_lon1D[ilatmid], dis_lat1D = pimplM->dis_lat1D;
	float grd_dis_lon = grd_lon * dis_lon1D,  grd_dis_lat = grd_lat * dis_lat1D;
	float grd_semidiag = sqrt( (grd_dis_lon * grd_dis_lon) + (grd_dis_lat * grd_dis_lat) );

	/* ellipse parameters */
	float Nmin = 3.; //(2 ~ 20?) Don't know much about sw kernel
	// define ellipse
	float f = dis*0.5; // known values
	float amax = f+lambda/(2.*Nmin);// asqrmax = amax*amax;
	//float bsqrmax = amax*amax - f*f; //maximum affective bsquare from Nmin
	//float bmax = sqrt(bsqrmax); // maximum affective width in the perpendicular direction
	float dab = amax - f, dab2 = dab*2.;
	float max_2a = 2. * (amax + grd_semidiag) + 10.;
	float max_esti = 2. * amax + 40.;


	float weit = 0., datasum = 0.;
	//float Nhaf = 12.;
	float alpha = -1.125 / (dab*dab); //- 0.5 / (dab*2.*0.33 * dab*2.*0.33);
	float dismax = 0.;
	auto fDist_ptr = acc ? &Path<float>::Dist : &Path<float>::DistF;
	for(int irow=0; irow<dataM.NumRows(); irow++) {
		for(int icol=0; icol<dataM.NumCols(); icol++) {
			// distances from (irow, icol) to src/rec
			float loncur = lonmin+irow*grd_lon, latcur = latmin+icol*grd_lat;
			float disEsrc = pimplM->estimate_dist( src, Point<float>(loncur,latcur) );
			float disErec = pimplM->estimate_dist( rec, Point<float>(loncur,latcur) );
			if( disEsrc + disErec > max_2a ) continue; // 20. for estimating error
			for( const auto& dpcur : dataM(irow, icol) ) {
				if( dpcur.Data() == NaN ) continue;
				//distance from dpcur to src/rec;
				float dis_src = dpcur.Dis(); //pimplM->estimate_dist(src, dpcur);
				float dis_rec = pimplM->estimate_dist(rec, dpcur);
				if( dis_src+dis_rec > max_esti ) continue; // 2.*dab == hdis * 3.
				dis_rec = (Path<float>(rec, dpcur).*fDist_ptr)();
				//calc_dist(src.Lat(), src.Lon(), dpcur.Lat(), dpcur.Lon(), &dis_src);
				//calc_dist(rec.Lat(), rec.Lon(), dpcur.Lat(), dpcur.Lon(), &dis_rec);
				float dis_ellip = dis_src + dis_rec - dis; // dis == 2.*f
				if( dis_ellip > dab2 ) continue; // 2.*dab == hdis * 3.
				if( dismax < dpcur.Dis() ) dismax = dpcur.Dis();
				float weight = exp( alpha * dis_ellip * dis_ellip );
				//std::cerr<<(Point<float>)dpcur<<" "<<weight<<"   "<<src<<"  "<<rec<<std::endl;
				weit += weight;
				datasum += (dpcur.Data() * weight);
			}
		}
	}
	if(weit==0.) {
		datasum = NaN; perc = 0.;
	} else {
		datasum /= weit;
		perc = dismax>dis ? 1 : dismax/dis;
	}
	//return datasum;
	return DataPoint<float>(rec, datasum, dis);

}


/* ------------ compute average along the path src-rec weighted by the reciprocal of the map value ------------ */
DataPoint<float> Map::PathAverage_Reci(Point<float> rec, float& perc, const float lambda, const bool acc) {
	// check source
	if( src == Point<float>() )
		throw ErrorM::BadParam(FuncName, "invalid src location");

	float dis = Path<float>(src, rec).Dist();
	if( lambda < 0. ) {
		throw ErrorM::BadParam(FuncName, "negative lambda");
	} else if ( lambda == 0. ) { // call direct ray tracing if lambda == 0
		//if( !pimplM->isReg )
		//	throw ErrorM::BadParam(FuncName, "lambda cannot be 0 when map is irregular");
		int N = 0; float avg = 0.;
		TraceGCP(src, rec, trace_step, [&](float val) { avg += 1./val; N++;} );
		avg = N / avg; perc = 1.;
		return DataPoint<float>(rec, avg, dis);
	}

	// references
	const auto& dataM = pimplM->dataM1;
	float lonmin = pimplM->lonmin, latmin = pimplM->latmin;
	float grd_lon = pimplM->grd1_lon, grd_lat = pimplM->grd1_lat;

	/* rec parameters */
	size_t idlmax = pimplM->dis_lon1D.size() - 1;
	int ilatmid = (int)floor( ( (rec.Lat()+src.Lat()) * 0.5 - latmin ) / grd_lat + 0.5 );
	if( ilatmid > idlmax ) ilatmid = idlmax;
	else if( ilatmid < 0 ) ilatmid = 0;
	float dis_lon1D = pimplM->dis_lon1D[ilatmid], dis_lat1D = pimplM->dis_lat1D;
	float grd_dis_lon = grd_lon * dis_lon1D,  grd_dis_lat = grd_lat * dis_lat1D;
	float grd_semidiag = sqrt( (grd_dis_lon * grd_dis_lon) + (grd_dis_lat * grd_dis_lat) );
	//std::cerr<<grd_dis_lon<<" "<<grd_dis_lat<<" "<<grd_semidiag<<"\n";

	/* ellipse parameters */
	float Nmin = 3.; //(2 ~ 20?) Don't know much about sw kernel; note that the data in the zone are weighted with a smaller hlaf length, so the effective N is larger!!!
	// define ellipse
	float f = dis*0.5; // known values
	float dab = lambda/(2.*Nmin), dab2 = dab*2.;
	float amax = f+dab;// asqrmax = amax*amax;
	//float bsqrmax = amax*amax - f*f; //maximum affective bsquare from Nmin
	//float bmax = sqrt(bsqrmax); // maximum affective width in the perpendicular direction
	float max_2a = 2. * (amax + grd_semidiag) + 10.;
	float max_esti = 2. * amax + 40.;


	float weit = 0., datasum = 0.;
	//float Nhaf = 12.;
	float alpha = -1.125 / (dab*dab); //- 0.5 / (dab*2.*0.33 * dab*2.*0.33);
	float dismax = 0., dismin = 99999.;
	//std::ofstream fout;
	//if( !outname.empty() ) fout.open(outname);
	auto fDist_ptr = acc ? &Path<float>::Dist : &Path<float>::DistF;
	//auto fDist_ptr = &Path<float>::Dist;
	for(int irow=0; irow<dataM.NumRows(); irow++) {
		for(int icol=0; icol<dataM.NumCols(); icol++) {
			// distances from (irow, icol) to src/rec
			float loncur = lonmin+irow*grd_lon, latcur = latmin+icol*grd_lat;
			float disEsrc = pimplM->estimate_dist( src, Point<float>(loncur,latcur) );
			float disErec = pimplM->estimate_dist( rec, Point<float>(loncur,latcur) );
			//float disEsrc1 = Path<float>( src, Point<float>(loncur,latcur) ).Dist();
			//float disErec1 = Path<float>( rec, Point<float>(loncur,latcur) ).Dist();
			//std::cerr<<Point<float>(loncur,latcur)<<" "<<rec<<"   "<<disEsrc<<" "<<disEsrc1<<"   "<<disErec<<" "<<disErec1<<"   "<<max_2a<<"   "<<grd_lon<<"\n";
//exit(0);
			if( disEsrc + disErec > max_2a ) continue; // 20. for estimating error
			for( const auto& dpcur : dataM(irow, icol) ) {
				if( dpcur.Data() == NaN ) continue;
				//distance from dpcur to src/rec;
				float dis_src = dpcur.Dis(); //pimplM->estimate_dist(src, dpcur);
				//float dis_rec1 = Path<float>(rec, dpcur).Dist();
				//std::cerr<<(Point<float>)src<<" "<<(Point<float>)dpcur<<"   "<<dis_src<<" "<<dis_rec<<" "<<dis_rec1<<"   "<<max_esti<<"   "<<grd_lon<<"\n";
				float dis_rec = pimplM->estimate_dist(rec, dpcur);
				if( dis_src+dis_rec > max_esti ) continue; // 2.*dab == hdis * 3.
				dis_rec = (Path<float>(rec, dpcur).*fDist_ptr)();
				//calc_dist(src.Lat(), src.Lon(), dpcur.Lat(), dpcur.Lon(), &dis_src);
				//calc_dist(rec.Lat(), rec.Lon(), dpcur.Lat(), dpcur.Lon(), &dis_rec);
				float dis_ellip = dis_src + dis_rec - dis; // dis == 2.*f
				if( dis_ellip > dab2 ) continue; // 2.*dab == hdis * 3.
				if( dismax < dpcur.Dis() ) dismax = dpcur.Dis();
				if( dismin > dpcur.Dis() ) dismin = dpcur.Dis();
				float weight = exp( alpha * dis_ellip * dis_ellip );
				if( weight < 0.01 ) continue;
				//std::cerr<<(Point<float>)dpcur<<" "<<weight<<"   "<<src<<"  "<<rec<<std::endl;
				weit += weight;
				//fout<<static_cast< Point<float> >(dpcur)<<" "<<dpcur.Data()<<" "<<weight<<std::endl;
				datasum += ( weight / dpcur.Data() );
				//std::cerr<<dpcur<<std::endl;
				//if( dpcur.Data() != dpcur.Data() ) std::cerr<<"nan!!! "<<dpcur<<"\n";
			}
		}
	}
	//std::cerr<<"datasum = "<<datasum<<"   weit = "<<weit<<"   dismax = "<<dismax<<std::endl;
	if(weit==0.) datasum = NaN;
	else datasum = weit/datasum;
	perc = dismax - dismin;
	perc = perc>=dis ? 1 : perc/dis;
	//return datasum;
	return DataPoint<float>(rec, datasum, dis);

}

