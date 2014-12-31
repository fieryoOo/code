// input:  map of points
// output: map of percentage on each grid

#include "Map.h"
#include <iostream>
#include <fstream>
#include <vector>

int main( int argc, char* argv[] ) {
	if( argc != 9 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [in_location_map (x y)] [xmin] [xmax] [gridlen_x] [ymin] [ymax] [gridlen_y] [outname]"<<std::endl;
		exit(-1);
	}

	// read in parameters
	float xmin = atof(argv[2]), xmax = atof(argv[3]), grdx = atof(argv[4]);
	float ymin = atof(argv[5]), ymax = atof(argv[6]), grdy = atof(argv[7]);
	if( grdx<=0 || grdy<=0 ) {
		std::cerr<<"Error(main): Invalid(non-positive) x/y grid length!"<<std::endl;
		exit(-2);
	}
	float hgrdx = grdx*0.5, hgrdy = grdy*0.5;

	// read in location map
	Map mapin(argv[1], hgrdx, hgrdy);

	// compute Number of points on each grid point
	std::vector< DataPoint<float> > Vout;
	//float lonmin = mapin.LonMin(), lonmax = mapin.LonMax();
	//float latmin = mapin.LatMin(), latmax = mapin.LatMax();
	for( float lon=xmin; lon<xmax; lon+=grdx )
		for( float lat=ymin; lat<ymax; lat+=grdy ) {
			float loneff, lateff;
			float Npts = mapin.NumberOfPoints( Point<float>(lon, lat), hgrdx, hgrdy, loneff, lateff );
			if( Npts==0 ) { loneff = lon; lateff = lat; }
			Vout.push_back( DataPoint<float>(loneff, lateff, Npts) );
		}

	// output
	std::ofstream fout(argv[8]);
	if( ! fout ) {
		std::cerr<<"Error(main): Cannot open file "<<argv[8]<<" to write!"<<std::endl;
		exit(-2);
	}
	for( const auto& dp : Vout ) fout<<dp<<"\n";
	fout.close();

	return 0;
}
