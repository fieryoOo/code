#ifndef BSPLINE_H
#define BSPLINE_H

#include "Curve.h"
#include "datatable.h"
#include "pspline.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
//#include <limits>


//template <class T>
class BSpline {
public:
	BSpline( const std::string& fname, const float penaltyin = 0., const int stepfactorin = 10 ) 
		: curve(fname), penalty(penaltyin), stepfactor(stepfactorin) {}

	BSpline( Curve<Point> curvein, const float penaltyin = 0., const int stepfactorin = 10 ) 
		: curve(curvein), penalty(penaltyin), stepfactor(stepfactorin) {}

	void Evaluate( Curve<Point>& curveout ) {
		using namespace MultivariateSplines;
		DenseVector loc(1);
      DataTable datain;
      Point dataP;
      for( curve.rewind(); curve.get(dataP); curve.next() ) {
			loc(0) = dataP.x;
			datain.addSample(loc, dataP.y);
		}

		// build P-spline (with penalty)
		PSpline pspline(datain, penalty);

		// output step
      float xmin = curve.xmin(), xmax = curve.xmax();
		float step = (xmax-xmin) / (datain.getNumSamples()*stepfactor);
		// output
		curveout.clear();
		curveout.reserve( (int)((xmax-xmin)/step) + 2 );
      for(float x=xmin; x<xmax; x+=step) {
         loc(0) = x;
			curveout.push_back( x, (float)pspline.eval(loc) );
      }

	}

/*
	void Load( const std::string& fname ) {
      // read from file
      std::ifstream fin( fname );
      if( ! fin )
			throw std::runtime_error("Error(BSplineLoad): cannot write to file " + fname);
	}
*/

private:
	float penalty = 0.;
	int stepfactor = 10;
	Curve<Point> curve;
};

#endif
