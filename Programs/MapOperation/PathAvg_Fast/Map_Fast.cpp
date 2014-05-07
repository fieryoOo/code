#include "Map_Fast.h"
#include "DisAzi.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <errno.h>

//int calc_dist(double lati1, double long1, double lati2, double long2, double *dist);

//int calc_azimuth(double lati1, double long1, double lati2, double long2, double *alpha1);

struct Map::Mimpl {
   std::string fname;
   Point<float> src;
   int ndis, ndisa;
   float grd_dis;
   float dis0, dismax;
   float disa0, disamax;
   std::vector< DataPoint<float> > dataV;
   Array2D< std::vector< DataPoint<float> > > dataM;


   Mimpl( const char* inname )
      : fname(inname) {
      dismax = disamax = -12345.;
   }

   Mimpl( const char* inname, const Point<float>& srcin )
      : fname(inname), src(srcin) {
      dismax = disamax = -12345.;
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
      if( autosrc ) {
	 std::string line;
	 if( ! std::getline(fin, line) ) return false;
	 float lon, lat, data;
	 sscanf(line.c_str(), "%f %f %f", &lon, &lat, &data);
	 lon = lon>90. ? lon-90. : lon+270.;
	 lat = lat>0 ? lat-90. : lat+90.;
	 src = Point<float>(lon, lat);
	 fin.seekg(0);
	 //dataV.push_back( DataPoint<float>(0., 0., data) );
      }
      for(std::string line; std::getline(fin, line); ) {
	 float lon, lat, data;
	 sscanf(line.c_str(), "%f %f %f", &lon, &lat, &data);
	 Path<float> pathcur(src, Point<float>(lon,lat));
	 //double dis, azi;
	 //calc_dist( src.Lat(), src.Lon(), lat, lon, &dis );
	 //calc_azimuth( src.Lat(), src.Lon(), lat, lon, &azi );
	 dataV.push_back( DataPoint<float>(pathcur.Dist(), pathcur.Azi1(), data) );
      }
      fin.close();
      if( dataV.size() == 0 ) return false;
      //std::cout<<"### "<<dataV.size()<<" map points loaded. ###"<<std::endl;
      /* map boundaries */
      for(size_t i=0; i<dataV.size(); i++) {
	 if( dismax < dataV[i].Dis() ) dismax = dataV[i].Dis();
	 if( disamax < dataV[i].Disa() ) disamax = dataV[i].Disa();
	 //std::cout<<dataV.at(i)<<std::endl;
      }
      return true;
   }

   void Hash() {
      if( dataV.size() == 0 ) return;
      grd_dis = 100.;
      dis0 = 0.;
      disa0 = -grd_dis * 2.;
      ndis = (int)ceil( (dismax-dis0) / grd_dis ) + 1;
      ndisa = (int)ceil( (disamax - disa0) / grd_dis ) + 1;
      dataM.Resize(ndis, ndisa);
      for(size_t i=0; i<dataV.size(); i++) {
	 int idis = (int)floor((dataV[i].Dis()-dis0) / grd_dis + 0.5);
	 int idisa = (int)floor((dataV[i].Disa()-disa0) / grd_dis + 0.5);
	 dataM(idis, idisa).push_back(dataV[i]);
	 // fill azimuth<0
	 DataPoint<float> dptmp(dataV[i].Dis(), dataV[i].Azi()-360., dataV[i].Data() );
	 idisa = (int)floor((dptmp.Disa()-disa0) / grd_dis + 0.5);
	 if( idisa >= 0 ) { dataM(idis, idisa).push_back(dptmp); }
	 // fill azimuth>360
	 dptmp = DataPoint<float>(dataV[i].Dis(), dataV[i].Azi()+360., dataV[i].Data() );
	 idisa = (int)floor((dptmp.Disa()-disa0) / grd_dis + 0.5);
	 if( idisa < ndisa ){ dataM(idis, idisa).push_back(dptmp); }
      }
/*
      // output
      for(int irow=0; irow<dataM.NumRows(); irow++) {
	 for(int icol=0; icol<dataM.NumCols(); icol++) {
	    for(size_t idata=0; idata<dataM(irow, icol).size(); idata++) {
		std::cerr<<dataM(irow, icol).at(idata)<<std::endl;;
	    }
	 }
      }
*/
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
   Point<float> src = pimplM->src;
   // compute dis and disa of Point rec
   Path<float> pathrec(src, rec);
   //double dis, azi;
   //calc_dist( src.Lat(), src.Lon(), rec.Lat(), rec.Lon(), &dis );
   //calc_azimuth( src.Lat(), src.Lon(), rec.Lat(), rec.Lon(), &azi );
   DataPoint<float> DPrec(pathrec.Dist(), pathrec.Azi1());

   Array2D< std::vector< DataPoint<float> > >& dataM = pimplM->dataM;
   float dis0 = pimplM->dis0, disa0 = pimplM->disa0, grd_dis = pimplM->grd_dis;

   // define computation area
   float dismax = hdis * 2.5, dismax_s = dismax*dismax;
   int rowmin = (int)floor((DPrec.Dis()-dismax-dis0) / grd_dis + 0.5);
   if( rowmin < 0 ) rowmin = 0;
   int rowmax = (int)floor((DPrec.Dis()+dismax-dis0) / grd_dis + 0.5) + 1;
   if( rowmax > dataM.NumRows() ) rowmax = dataM.NumRows();
   int colmin = (int)floor((DPrec.Disa()-dismax-disa0) / grd_dis + 0.5);
   if( colmin < 0 ) colmin = 0;
   int colmax = (int)floor((DPrec.Disa()+dismax-disa0) / grd_dis + 0.5) + 1;
   if( colmax > dataM.NumCols() ) colmax = dataM.NumCols();

   weit = 0.;
   float alpha = -0.5 / (hdis*hdis), datasum = 0.;
   for(int irow=rowmin; irow<rowmax; irow++) {
      for(int icol=colmin; icol<colmax; icol++) {
         for(size_t idata=0; idata<dataM(irow, icol).size(); idata++) {
	    DataPoint<float> dpcur = dataM(irow, icol)[idata];
	    //distance from dataM(irow, icol).at(idata) to DPrec.
	    float xdis = DPrec.Dis() - dpcur.Dis(), ydis = DPrec.Disa() - dpcur.Disa();
	    float disc_s = xdis*xdis + ydis*ydis;
//std::cerr<<dpcur.Data()<<"  "<<disc_s<<"  "<<dpcur.Azi()<<" "<<dpcur.Dis()<<"  "<<DPrec.Azi()<<" "<<DPrec.Dis()<<std::endl;
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
   Point<float> src = pimplM->src;
   // compute dis and disa of Point rec
   //double dis, azi;
   //calc_dist( src.Lat(), src.Lon(), rec.Lat(), rec.Lon(), &dis );
   //calc_azimuth( src.Lat(), src.Lon(), rec.Lat(), rec.Lon(), &azi );
   Path<float> pathrec(src, rec);
   double dis = pathrec.Dist();
   DataPoint<float> DPrec(dis, pathrec.Azi1());

   float Nmin = 3.; //(2 ~ 20?) Don't know much about sw kernel
   // define ellipse
   float f = dis*0.5; // known values
   float amax = f+lamda/(2.*Nmin), asqrmax = amax*amax;
   float bsqrmax = amax*amax - f*f; //maximum affective bsquare from Nmin
   float bmax = sqrt(bsqrmax); // maximum affective width in the perpendicular direction
   float dab = amax - f;

   Array2D< std::vector< DataPoint<float> > >& dataM = pimplM->dataM;
   float dis0 = pimplM->dis0, disa0 = pimplM->disa0, grd_dis = pimplM->grd_dis;

   int rowmax = (int)floor((DPrec.Dis()+dab-dis0) / grd_dis + 0.5) + 1;
   if( rowmax > dataM.NumRows() ) rowmax = dataM.NumRows();
   int colmax = (int)floor( ( std::max(DPrec.Disa(),(float)0.)+bmax-disa0 ) / grd_dis + 0.5 ) + 1;
   if( colmax > dataM.NumCols() ) colmax = dataM.NumCols();
   int colmin = (int)floor( ( std::min(DPrec.Disa(),(float)0.)-bmax-disa0 ) / grd_dis + 0.5) - 1;
   if( colmin < 0 ) colmin = 0;
   float weit = 0., datasum = 0.;
   //float Nhaf = 12.;
   float dismax = 0.;
   for(int irow=0; irow<rowmax; irow++) {
      for(int icol=colmin; icol<colmax; icol++) {
	 // distance from (irow, icol) to the line(src, DPrec)
	 float disP = (disa0 + grd_dis * icol) - (DPrec.Disa() * (dis0+grd_dis*irow)/DPrec.Dis());
	 if( fabs(disP) > bmax+grd_dis ) continue;
//std::cerr<<disa0+grd_dis*icol<<" "<<dis0+grd_dis*irow<<std::endl;
         for(size_t idata=0; idata<dataM(irow, icol).size(); idata++) {
	    DataPoint<float> dpcur = dataM(irow, icol)[idata];
	    //distance from dataM(irow, icol).at(idata) to line(src, DPrec);
	    float xdat = dpcur.Dis() - f;
	    if( fabs(xdat) > amax ) continue;
	    float ydat = dpcur.Disa() - (DPrec.Disa() * dpcur.Dis()/DPrec.Dis());
	    float ymax = sqrt( 1-(xdat*xdat/asqrmax) ) * bmax;
	    if( fabs(ydat) > ymax ) continue;
	    if( dismax < dpcur.Dis() ) dismax = dpcur.Dis();
	    float factor = ydat / (ymax * 0.33);
	    float weight = exp(-0.5 * factor * factor);
	    //std::cerr<<dpcur.Dis()<<" "<<dpcur.Disa()<<" "<<weight<<"   "<<src<<"  "<<rec<<std::endl;
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
   DPrec.Data() = datasum;
   return DPrec;
 
}

