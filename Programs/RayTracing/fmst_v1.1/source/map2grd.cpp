#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

class PointData {
public:
   float lon, lat;
   float dat;

   /* con/destructors and operators */
   PointData() {}
   PointData( const char *info ) { Load(info); }
   PointData( float lonin, float latin, float datin ) 
      : lon(lonin), lat(latin), dat(datin) {}
   ~PointData() {}
   friend std::ostream& operator<< (std::ostream& o, const PointData& pd) { 
      o<<"( "<<pd.lon<<", "<<pd.lat<<" ): "<<pd.dat; 
      return o; 
   }

   /* load point from string */
   void Load( const char *info ) { 
      sscanf( info, "%f %f %f", &lon, &lat, &dat ); 
      if( lon < 0 ) lon += 360.;
   }

   /* location comparisons */
   bool IsWestOf( PointData pd2 ) { return (lon < pd2.lon); }
   bool IsNorthOf( PointData pd2 ) { return (lat > pd2.lat); }

};

int main( int argc, char *argv[] ) {
   /* input */
   if( argc != 5 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [input map] [lon step] [lat step] [output file]"<<std::endl;
      exit(-1);
   }
   
   /* open input map */
   std::ifstream fin(argv[1]);
   if( ! fin ) {
      std::cerr<<"Cannot read from file "<<argv[1]<<std::endl;
      exit(0);
   }

   /* read in point data */
   std::vector<PointData> Map;
   for( std::string line; std::getline(fin, line); ){
      PointData pdtmp( line.c_str() );
      Map.push_back(pdtmp);
   }  
   fin.close();
   if( Map.size() <= 0 ) {
      std::cerr<<"Empty map!"<<std::endl;
      exit(0);
   }

   /* search for boundaries */
   float lonstep = atof(argv[2]), latstep = atof(argv[3]);
   if( lonstep<=0 || latstep<=0 ) {
      std::cerr<<"positive number expected for lon/lat steps!"<<std::endl;
      exit(0);
   }
   PointData UL( Map.at(0) ); // upper-left corner
   PointData LR( Map.at(0) ); // lower-right corner
   for( int i=0; i<Map.size(); i++ ) {
      PointData &pdcurrent = Map.at(i);
      if( pdcurrent.IsWestOf(UL) ) UL.lon = pdcurrent.lon;
      else if( LR.IsWestOf(pdcurrent) ) LR.lon = pdcurrent.lon;
      if( pdcurrent.IsNorthOf(UL) ) UL.lat = pdcurrent.lat;
      else if( LR.IsNorthOf(pdcurrent) ) LR.lat = pdcurrent.lat;
   }
   std::cout<<"Upper-left: "<<UL<<std::endl;
   std::cout<<"Lower-right: "<<LR<<std::endl;

   /* compute data average */
   float mean = 0.;
   for( int i=0; i<Map.size(); i++ ) mean += Map.at(i).dat;
   mean /= Map.size();
   std::cout<<"data mean = "<<mean<<std::endl;

   /* data matrix */
   int npts_lon = (int)floor( (LR.lon - UL.lon)/lonstep + 0.5 ) + 1;
   int npts_lat = (int)floor( (UL.lat - LR.lat)/latstep + 0.5 ) + 1;
   std::vector<float> DataM;
   DataM.resize( npts_lon * npts_lat );
   std::fill( DataM.begin(), DataM.end(), mean);
   for( int i=0; i<Map.size(); i++ ) {
      PointData &pdcurrent = Map.at(i);
      int idxlon = (int)floor( (pdcurrent.lon - UL.lon) / lonstep + 0.5 );
      int idxlat = (int)floor( (UL.lat - pdcurrent.lat) / latstep + 0.5 );
      DataM.at( idxlon * npts_lat + idxlat ) = pdcurrent.dat;
      //std::cerr<<idxlat<<" "<<idxlon<<" "<<pdcurrent<<std::endl;
   }

   /* output in the format required by fm2dss */
   std::ofstream fout(argv[4]);
   if( ! fout ) {
      std::cerr<<"Cannot write to file "<<argv[3]<<std::endl;
      exit(0);
   }
   fout << npts_lat-2 << " " << npts_lon-2 << std::endl; // leave the outer boundary as 'cushion nodes'
   fout << UL.lat-latstep << " " << UL.lon+lonstep << std::endl; // computational grid starts from the second point on each direction
   fout << latstep << " " << lonstep << std::endl; 
   for( int i=0; i<DataM.size(); i++ ) fout << DataM.at(i) << std::endl;
   fout.close();

   return 0;
}

