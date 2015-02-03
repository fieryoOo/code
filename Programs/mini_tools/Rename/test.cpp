#include "MyOMP.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

int main(int argc, char* argv[]) {
	// test arguments
	if( argc != 4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [number1] [number2] [number3]"<<std::endl;
		return -1;
	}
	float f1 = atof(argv[1]), f2 = atof(argv[2]), f3 = atof(argv[3]);
	std::cout<<"Input arguments are "<<f1<<" "<<f2<<" "<<f3<<std::endl;

	// test stdin
	std::string linein;
	std::getline(std::cin, linein);
	std::cout<<"stdin is "<<linein<<std::endl;

	// test openmp
	const int nraw = 300, ncol = 300;
	float mat[nraw][ncol];
	int Ghalfdis = 80;
	float alpha = -0.5 / (Ghalfdis*Ghalfdis);
	#pragma omp parallel for
	for(int iraw=0; iraw<nraw; iraw++) {
		for(int icol=0; icol<ncol; icol++) {
			int disx = icol - 100;
			int disy = iraw - 200;
			int disS = disx*disx + disy*disy;
			mat[iraw][icol] = exp(alpha * disS);
		}
	}
	mat[150][280] = -1.;

	return 0;
}
