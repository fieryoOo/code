#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

struct Data2C {
	float x = NaN, y = NaN;

	Data2C( const std::string& input ) {
		if( sscanf(input.c_str(), "%f %f", &x, &y) != 2 )
			x = y = NaN;
	}

	friend std::ostream& operator<< ( std::ostream& o, const Data2C& data ) {
		o << data.x << "\t" << data.y;
		return o;
	}

	bool isValid() {
		return (x!=NaN && y!=NaN);
	}

protected:
	static constexpr float NaN = -12345.;
};

inline bool CompareABS( const Data2C& d1, const Data2C& d2 ) {
	return ( fabs(d1.y) < fabs(d2.y) );
}


int main( int argc, char* argv[] ) {
	// check inputs
	if( argc != 4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [in_file] [max amp] [out_file]"<<std::endl;
		exit(-1);
	}

	// read in data file
	std::ifstream fin(argv[1]);
	if( ! fin ) {
		std::cerr<<"Error(main): cannot read from file "<<argv[1]<<std::endl;
		exit(-2);
	}
	std::vector<Data2C> dataV;
	for(std::string line; std::getline(fin, line); ) {
		Data2C data_cur(line);
		if( data_cur.isValid() ) // check if read is succeed
			dataV.push_back( data_cur );
	}
	fin.close();
	// empty data?
	if( dataV.size() == 0 ) {
		std::cerr<<"Error(main): empty data!"<<std::endl;
		exit(-3);
	} else {
		std::cout<<dataV.size()<<" data points read in"<<std::endl;
	}

	// find max absolute
	auto iter_max = std::max_element( dataV.begin(), dataV.end(), CompareABS );
	float norm_factor = fabs( (*iter_max).y );
	if( norm_factor == 0. ) {
		std::cerr<<"Error(main): all zero inputs!"<<std::endl;
		exit(-3);
	}
	norm_factor = atof(argv[2]) / norm_factor;

	// normalize
	for( auto& data_cur : dataV )
		data_cur.y *= norm_factor;

	// output
	std::ofstream fout(argv[3]);
	if( ! fout ) {
		std::cerr<<"Error(main): cannot write to file "<<argv[3]<<std::endl;
		exit(-2);
	}
	for( const auto& data_cur : dataV )
		fout << data_cur << "\n";
	
	
	return 0;
}
