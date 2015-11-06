#include "spline.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include <stdexcept>

int main( int argc, char* argv[] ) {
	if( argc != 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [input file (x, y)]"<<std::endl;
		return -1;
	}
	
	// open/check file
	std::string fname(argv[1]);
	std::ifstream fin(fname);
	std::string funcname(__FUNCTION__);
	if( ! fin )
		throw std::runtime_error("Error("+funcname+"): IO failed on "+fname);
	// load data and record x range
	std::vector<double> xV, yV;
	float xmin = std::numeric_limits<float>::max();
	float xmax = std::numeric_limits<float>::lowest();
	for(std::string line; std::getline(fin, line);) {
		double x, y;
		if( sscanf(line.c_str(), "%lf %lf", &x, &y) != 2 ) {
			std::cerr<<"Warning("<<funcname<<"): input format error in file "<<fname<<std::endl;
			continue;
		}
		xV.push_back(x);
		yV.push_back(y);
		if( x < xmin ) {
			xmin = x;
		} else if( x > xmax ) {
			xmax = x;
		}
	}
	// compute spline
	tk::spline s; s.set_points(xV, yV);
	// output
	std::ofstream fout(fname+"_spline");
	int nstep = 100;
	float step = (xmax - xmin)/(nstep-1);
	for(int i=-1; i<nstep+1; i++) {
		double x = xmin + i*step;
		fout<<x<<" "<<s(x)<<"\n";
	}
	return 0;
}
