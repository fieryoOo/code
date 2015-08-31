#include "BSpline.h"
#include "datatable.h"
#include "pspline.h"
//#include "bspline.h"
//#include "rbfspline.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>

//using namespace MultivariateSplines;

int main( int argc, char* argv[] ) {
	if( argc != 3 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [input file (x, y)] [smooth penalty(0. - 1.?)]"<<std::endl;
		return -1;
	}

	try {
		BSpline bspline( argv[1], atof(argv[2]) );
		Curve<PointC> curveout;
		bspline.Evaluate( curveout );

		// output
		std::string outname(argv[1]);
		outname += "_spline";
		curveout.Write(outname);

	} catch ( MultivariateSplines::Exception& e ) {
		std::cerr<<"MSException(main): "<<e.what()<<std::endl;
	} catch ( std::exception& e ) {
		std::cerr<<"stdException(main): "<<e.what()<<std::endl;
	} catch (...) {
		std::cerr<<"Error(main): unknown exception detected!"<<std::endl;
	}


	return 0;
}
