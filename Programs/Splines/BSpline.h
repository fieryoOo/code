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

	BSpline( Curve<PointC> curvein, const float penaltyin = 0., const int stepfactorin = 10 ) 
		: curve(curvein), penalty(penaltyin), stepfactor(stepfactorin) {}

	void Evaluate( Curve<PointC>& curveout ) {
		using namespace MultivariateSplines;
		DenseVector loc(1);
      DataTable datain;
      PointC dataP;
      for( curve.rewind(); curve.get(dataP); curve.next() ) {
			loc(0) = dataP.x;
			datain.addSample(loc, dataP.y);
		}

		// build P-spline (with penalty)
		PSpline pspline(datain, penalty);

		// output step
      double xmin = curve.xmin(), xmax = curve.xmax();
		int nstep = datain.getNumSamples()*stepfactor;
		double step = (xmax-xmin) / nstep;
		xmin += step; xmax -= step;	// prevent out-of-bound exception
		step = (xmax-xmin) / nstep;
		
		// output
		curveout.clear();
		curveout.reserve( (int)((xmax-xmin)/step) + 2 );
      for(double x=xmin; x<xmax; x+=step) {
         loc(0) = x;
			//std::cerr<<"Evaluate at "<<x<<std::endl;
			curveout.push_back( x, (double)pspline.eval(loc) );
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
	Curve<PointC> curve;
};

#endif
