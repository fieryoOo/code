#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>


bool ExcludeBad( std::vector<float>& data, std::vector<float>& weight, std::vector<std::string>& fileV, float stdfactor=2. ) {
	if( data.size() != weight.size() ) return false;
	// compute mean 1
	float V1 = 0., mean1 = 0.;
	for(int i=0; i<data.size(); i++) {
		mean1 += data[i] * weight[i];
		V1 += weight[i];
	}
	mean1 /= V1;
	// compute std1
	float V2 = 0., std1 = 0.;
	for(int i=0; i<data.size(); i++) {
		float ftmp = data[i]-mean1;
		std1 += ftmp * ftmp * weight[i];
		V2 += weight[i] * weight[i];
	}
	std1 = sqrt( std1 * V1 / (V1*V1-V2) );
	std::cerr<<mean1<<" "<<std1<<std::endl;
	// exclude larger-than-2.0sigma data
	float ubound = mean1 + stdfactor*std1, lbound = mean1 - stdfactor*std1;
	int i=0;
	while( i<data.size() ) {
		if( data[i]>ubound || data[i]<lbound ) {
			data.erase( data.begin() + i );
			weight.erase( weight.begin() + i );
			fileV.erase( fileV.begin() + i );
		} else {
			i++;
		}
	}

	return true;
}

int main( int argc, char* argv[] ) {
	/* check input */
	if( argc != 5 ) {
		std::cerr<<"Usage "<<argv[0]<<" [input_file] [data column] [std factor( exclude +/- factor*std )] [#iteration]"<<std::endl;
		exit(-1);
	}
	int icol=atoi(argv[2]);
	if( icol<=0 ) {
		std::cerr<<"Error(main): positive integer expected for data column number! icol = "<<icol<<std::endl;
		exit(0);
	}

	/* load in input file */
	std::ifstream fin(argv[1]);
	if( !fin ) {
		std::cerr<<"Error(main): Cannot read from file "<<argv[1]<<std::endl;
		exit(0);
	}
	std::vector<std::string> fileV;
	std::vector<float> dataV;
	for(std::string line; std::getline(fin, line); ) {
		std::string datastr;
		std::stringstream sin(line);
		for(int i=0; i<icol; i++) {
			//if(! std::getline(sin, datastr, ' ') ) {
			if(! (sin >> datastr) ) {
				std::cerr<<"Error(main): Invalid column number!"<<std::endl;
				exit(0);
			}
		}
		fileV.push_back(line);
		dataV.push_back(std::stof(datastr));
	}
	fin.close();

	/* exclude */
	int niter=atoi(argv[4]);
	if( niter<=0 ) {
		std::cerr<<"Error(main): positive integer expected for number of iterations! niter = "<<niter<<std::endl;
		exit(0);
	}
	float stdfactor=atof(argv[3]);
	std::vector<float> weit(dataV.size(), 1.);
	for(int i=0; i<niter; i++) ExcludeBad(dataV, weit, fileV, stdfactor);

	/* output */
	for(int i=0; i<fileV.size(); i++) std::cout<<fileV[i]<<std::endl;

	return 0;
}
