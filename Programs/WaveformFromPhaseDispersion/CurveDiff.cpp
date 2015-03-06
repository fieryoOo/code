#include "SacRec.h"
#include "Curve.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>


int main(int argc, char* argv[]) {
	if( argc!=4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [fcurve1 (x y)] [fcurve2] [fcurveout]"<<std::endl;
		exit(-1);
	}

	try {
		Curve<Point> c1(argv[1]), c2(argv[2]);
		Curve<Point> cout = c2 - c1;
		cout.Write(argv[3]);
		/*
		std::ofstream fout(argv[3]);
		if( ! fout )
			throw std::runtime_error(std::string("Error(main): cannot write to file ") + argv[3]);
		for( auto p2 : c2.dataV ) {
			float y1 = c1.Val(p2.per);
			p2.vel -= y1;
			fout<<p2<<"\n";
		}
		*/
	} catch (...) {
		std::cerr<<"Exception detected!\n";
		exit(-3);
	}
}

