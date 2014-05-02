#include "Map.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <errno.h>

int calc_dist(double lati1, double long1, double lati2, double long2, double *dist);

int calc_azimuth(double lati1, double long1, double lati2, double long2, double *alpha1);

struct Map::Mimpl {
   std::string fname;
   Point<float> src;
   int nlon, nlat;
   float grd_lon, grd_lat;
   float lonmin, lonmax;
   float latmin, latmax;

   double dis_lat1D;
   std::vector<double> dis_lon1D;

   std::vector< DataPoint<float> > dataV;
   Array2D< std::vector< DataPoint<float> > > dataM;


   Mimpl( const char* inname, const Point<float> srcin = Point<float>() )
      : fname(inname), src(srcin)
      , nlon(-12345), nlat(-12345), grd_lon(-12345.), grd_lat(-12345.) 
      , lonmin(-12345.), lonmax(-12345.), latmin(-12345.), latmax(-12345) {
   }


   bool Load( bool autosrc = false ) {
      /* check src location */
      if( ! autosrc ) {
	 if( src.Lon() == -12345. || src.Lat() == -12345. ) return false;
      }
      /* read in the map */
      std::ifstream fin(fname.c_str());
      if( ! fin ) return false;
      dataV.clear();
      /* read the 1st line */
      std::string line;
      if( ! std::getline(fin, line) ) return false;
      float lon, lat, data;
      sscanf(line.c_str(), "%f %f %f", &lon, &lat, &data);
      lonmin = lonmax = lon;
      latmin = latmax = lat;
      if( autosrc ) {
	 lon = lon>90. ? lon-90. : lon+270.;
	 lat = lat>0 ? lat-90. : lat+90.;
	 src = Point<float>(lon, lat);
      }
      fin.seekg(0); // rewind
      /* load in all data */
      for(std::string line; std::getline(fin, line); ) {
	 float lon, lat, data;
	 sscanf(line.c_str(), "%f %f %f", &lon, &lat, &data);
	 double dis;//, azi;
	 calc_dist( src.Lat(), src.Lon(), lat, lon, &dis );
	 //calc_azimuth( src.Lat(), src.Lon(), lat, lon, &azi );
	 dataV.push_back( DataPoint<float>(lon, lat, data, dis) );
      }
      fin.close();
      if( dataV.size() == 0 ) return false;
      //std::cout<<"### "<<dataV.size()<<" map points loaded. ###"<<std::endl;
      /* map boundaries */
      for(size_t i=0; i<dataV.size(); i++) {
	 DataPoint<float>& dpcur = dataV[i];
	 if( lonmax < dpcur.Lon() ) {
	    lonmax = dpcur.Lon();
	 } 
	 else if( lonmin > dpcur.Lon() ) {
	    lonmin = dpcur.Lon();
	 }
	 if( latmax < dpcur.Lat() ) {
	    latmax = dpcur.Lat();
	 } 
	 else if( latmin > dpcur.Lat() ) {
	    latmin = dpcur.Lat();
	 }
      }
      return true;
   }

   void Hash() {
      if( dataV.size() == 0 ) return;
      /* grid/hash size */
      //grd_lon = grd_lat = sqrt( (lonmax-lonmin) * (latmax-latmin) / sqrt((float)dataV.size()) );
      grd_lon = grd_lat = 1.;
      nlon = (int)ceil( (lonmax-lonmin) / grd_lon ) + 1;
      nlat = (int)ceil( (latmax-latmin) / grd_lat ) + 1;
      dataM.Resize(nlon, nlat);
      /* hash */
      for(size_t i=0; i<dataV.size(); i++) {
	 DataPoint<float>& dpcur = dataV[i];
	 int ilon = (int)floor((dpcur.Lon()-lonmin) / grd_lon + 0.5);
	 int ilat = (int)floor((dpcur.Lat()-latmin) / grd_lat + 0.5);
	 dataM(ilon, ilat).push_back(dpcur);
      }
      /* compute distance of 1 degree in lon/lat */
      dis_lat1D = 111.;
      dis_lon1D.resize( dataM.NumCols() );
      for(int ilat=0; ilat<dataM.NumCols(); ilat++) {
	 float latcur = latmin + ilat*grd_lat;
	 calc_dist( latcur, 0., latcur, 1., &(dis_lon1D[ilat]) );
      }
   }

   inline float estimate_dist(const Point<float>& p1, const Point<float>& p2) {
      float lon1 = p1.Lon(), lat1 = p1.Lat();
      if( lon1 < 0. ) lon1 += 360.;
      float lon2 = p2.Lon(), lat2 = p2.Lat();
      if( lon2 < 0. ) lon2 += 360.;

      int ilatmid = (int)floor( ( (lat1+lat2) * 0.5 - latmin ) / grd_lat + 0.5 );
      if( ilatmid > dataM.NumCols() ) ilatmid = dataM.NumCols();
      else if( ilatmid < 0 ) ilatmid = 0;

      float dis_lon = (lon1-lon2) * dis_lon1D[ilatmid];
      float dis_lat = (lat1-lat2) * dis_lat1D;
      return sqrt(dis_lon*dis_lon + dis_lat*dis_lat);
   }
};


/* -------------- con/destructors ----------------- */
Map::Map( const char *inname ) 
   : pimplM( new Mimpl(inname) ) {
   if( ! pimplM->Load( true ) ) { // Load with auto-source
      std::cerr<<"Error(Map::Map): Load failed!"<<std::endl;
      exit(0);
   }
   pimplM->Hash();
}

Map::Map( const char *inname, const Point<float>& srcin ) 
   : pimplM( new Mimpl(inname, srcin) ) {
   if( ! pimplM->Load( false ) ) { // Load without auto-source
      std::cerr<<"Error(Map::Map): Load failed!"<<std::endl;
      exit(0);
   }
   pimplM->Hash();
}


Map::~Map() {}



/* ------------ compute average value on the point rec ------------ */
float Map::PointAverage(Point<float> rec, float hdis, float& weit) {
   /* references */
   Array2D< std::vector< DataPoint<float> > >& dataM = pimplM->dataM;
   float lonmin = pimplM->lonmin, latmin = pimplM->latmin;
   float grdlon = pimplM->grd_lon, grdlat = pimplM->grd_lat;;

   // define computation area
   float dismax = hdis * 2.5, dismax_s = dismax*dismax;
   int ilatrec = (int)floor((rec.Lat()-latmin) / grdlat + 0.5);
   if( ilatrec > dataM.NumCols() ) ilatrec = dataM.NumCols();
   else if( ilatrec < 0 ) ilatrec = 0;
   double dis_lon1D = pimplM->dis_lon1D[ilatrec], dis_lat1D = pimplM->dis_lat1D;
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
         for(size_t idata=0; idata<dataM(irow, icol).size(); idata++) {
	    DataPoint<float> dpcur = dataM(irow, icol)[idata];
	    //distance from dataM(irow, icol).at(idata) to DPrec.
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
   if(weit==0.) datasum = -12345.;
   else datasum /= weit;

   return datasum;
 
}



/* ------------ compute average value along the path src-rec ------------ */
DataPoint<float> Map::PathAverage(Point<float> rec, float lamda, float& perc) {
   // references
   Array2D< std::vector< DataPoint<float> > >& dataM = pimplM->dataM;
   float lonmin = pimplM->lonmin, latmin = pimplM->latmin;
   float grd_lon = pimplM->grd_lon, grd_lat = pimplM->grd_lat;

   Point<float> src = pimplM->src;
   /* rec parameters */
   double dis;//, azi;
   calc_dist( src.Lat(), src.Lon(), rec.Lat(), rec.Lon(), &dis );
   //calc_azimuth( src.Lat(), src.Lon(), rec.Lat(), rec.Lon(), &azi );

   int ilatmid = (int)floor( ( (rec.Lat()+src.Lat()) * 0.5 - latmin ) / grd_lat + 0.5 );
   if( ilatmid > dataM.NumCols() ) ilatmid = dataM.NumCols();
   else if( ilatmid < 0 ) ilatmid = 0;
   double dis_lon1D = pimplM->dis_lon1D[ilatmid], dis_lat1D = pimplM->dis_lat1D;
   float grd_dis_lon = grd_lon * dis_lon1D,  grd_dis_lat = grd_lat * dis_lat1D;
   float grd_semidiag = sqrt( (grd_dis_lon * grd_dis_lon) + (grd_dis_lat * grd_dis_lat) );

   /* ellipse parameters */
   float Nmin = 3.; //(2 ~ 20?) Don't know much about sw kernel
   // define ellipse
   float f = dis*0.5; // known values
   float amax = f+lamda/(2.*Nmin);// asqrmax = amax*amax;
   //float bsqrmax = amax*amax - f*f; //maximum affective bsquare from Nmin
   //float bmax = sqrt(bsqrmax); // maximum affective width in the perpendicular direction
   float dab = amax - f, dab2 = dab*2.;
   float max_2a = 2. * (amax + grd_semidiag);
   float max_esti = 2.*amax+20.;


   float weit = 0., datasum = 0.;
   //float Nhaf = 12.;
   float alpha = -1.125 / (dab*dab); //- 0.5 / (dab*2.*0.33 * dab*2.*0.33);
   float dismax = 0.;
   for(int irow=0; irow<dataM.NumRows(); irow++) {
      for(int icol=0; icol<dataM.NumCols(); icol++) {
	 // distances from (irow, icol) to src/rec
	 float loncur = lonmin+irow*grd_lon, latcur = latmin+icol*grd_lat;
	 float disEsrc = pimplM->estimate_dist( src, Point<float>(loncur,latcur) );
	 float disErec = pimplM->estimate_dist( rec, Point<float>(loncur,latcur) );
	 if( disEsrc + disErec > max_2a+20. ) continue; // 20. for estimating error
         for(size_t idata=0; idata<dataM(irow, icol).size(); idata++) {
	    DataPoint<float> dpcur = dataM(irow, icol)[idata];
	    //distance from dataM(irow, icol).at(idata) to src/rec;
	    double dis_src = dpcur.Dis(); //pimplM->estimate_dist(src, dpcur);
	    double dis_rec = pimplM->estimate_dist(rec, dpcur);
	    if( dis_src+dis_rec > max_esti ) continue; // 2.*dab == hdis * 3.
	    //calc_dist(src.Lat(), src.Lon(), dpcur.Lat(), dpcur.Lon(), &dis_src);
	    calc_dist(rec.Lat(), rec.Lon(), dpcur.Lat(), dpcur.Lon(), &dis_rec);
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
   if(weit==0.) datasum = -12345.;
   else {
      datasum /= weit;
      perc = dismax>dis ? 1 : dismax/dis;
   }
   //return datasum;
   return DataPoint<float>(rec, datasum, dis);
 
}


/* ------------ compute average along the path src-rec weighted by the reciprocal of the map value ------------ */
DataPoint<float> Map::PathAverage_Reci(Point<float> rec, float lamda, float& perc) {
   // references
   Array2D< std::vector< DataPoint<float> > >& dataM = pimplM->dataM;
   float lonmin = pimplM->lonmin, latmin = pimplM->latmin;
   float grd_lon = pimplM->grd_lon, grd_lat = pimplM->grd_lat;

   Point<float> src = pimplM->src;
   /* rec parameters */
   double dis;//, azi;
   calc_dist( src.Lat(), src.Lon(), rec.Lat(), rec.Lon(), &dis );
   //calc_azimuth( src.Lat(), src.Lon(), rec.Lat(), rec.Lon(), &azi );

   int ilatmid = (int)floor( ( (rec.Lat()+src.Lat()) * 0.5 - latmin ) / grd_lat + 0.5 );
   if( ilatmid > dataM.NumCols() ) ilatmid = dataM.NumCols();
   else if( ilatmid < 0 ) ilatmid = 0;
   double dis_lon1D = pimplM->dis_lon1D[ilatmid], dis_lat1D = pimplM->dis_lat1D;
   float grd_dis_lon = grd_lon * dis_lon1D,  grd_dis_lat = grd_lat * dis_lat1D;
   float grd_semidiag = sqrt( (grd_dis_lon * grd_dis_lon) + (grd_dis_lat * grd_dis_lat) );

   /* ellipse parameters */
   float Nmin = 3.; //(2 ~ 20?) Don't know much about sw kernel
   // define ellipse
   float f = dis*0.5; // known values
   float amax = f+lamda/(2.*Nmin);// asqrmax = amax*amax;
   //float bsqrmax = amax*amax - f*f; //maximum affective bsquare from Nmin
   //float bmax = sqrt(bsqrmax); // maximum affective width in the perpendicular direction
   float dab = amax - f, dab2 = dab*2.;
   float max_2a = 2. * (amax + grd_semidiag);
   float max_esti = 2.*amax+20.;


   float weit = 0., datasum = 0.;
   //float Nhaf = 12.;
   float alpha = -1.125 / (dab*dab); //- 0.5 / (dab*2.*0.33 * dab*2.*0.33);
   float dismax = 0.;
   for(int irow=0; irow<dataM.NumRows(); irow++) {
      for(int icol=0; icol<dataM.NumCols(); icol++) {
	 // distances from (irow, icol) to src/rec
	 float loncur = lonmin+irow*grd_lon, latcur = latmin+icol*grd_lat;
	 float disEsrc = pimplM->estimate_dist( src, Point<float>(loncur,latcur) );
	 float disErec = pimplM->estimate_dist( rec, Point<float>(loncur,latcur) );
	 if( disEsrc + disErec > max_2a+20. ) continue; // 20. for estimating error
         for(size_t idata=0; idata<dataM(irow, icol).size(); idata++) {
	    DataPoint<float> dpcur = dataM(irow, icol)[idata];
	    //distance from dataM(irow, icol).at(idata) to src/rec;
	    double dis_src = dpcur.Dis(); //pimplM->estimate_dist(src, dpcur);
	    double dis_rec = pimplM->estimate_dist(rec, dpcur);
	    if( dis_src+dis_rec > max_esti ) continue; // 2.*dab == hdis * 3.
	    //calc_dist(src.Lat(), src.Lon(), dpcur.Lat(), dpcur.Lon(), &dis_src);
	    calc_dist(rec.Lat(), rec.Lon(), dpcur.Lat(), dpcur.Lon(), &dis_rec);
	    float dis_ellip = dis_src + dis_rec - dis; // dis == 2.*f
	    if( dis_ellip > dab2 ) continue; // 2.*dab == hdis * 3.
	    if( dismax < dpcur.Dis() ) dismax = dpcur.Dis();
	    float weight = exp( alpha * dis_ellip * dis_ellip );
	    //std::cerr<<(Point<float>)dpcur<<" "<<weight<<"   "<<src<<"  "<<rec<<std::endl;
	    weit += weight;
	    datasum += ( weight / dpcur.Data() );
         }
      }
   }
   if(weit==0.) datasum = -12345.;
   else datasum = weit/datasum;
   perc = dismax>dis ? 1 : dismax/dis;
   //return datasum;
   return DataPoint<float>(rec, datasum, dis);
 
}

