#include <SacRec.h>
#include <iostream>
#include <fstream>
#include <string>

int main( int argc, char* argv[] ) {
   if( argc != 3 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [sac_list] [sac_out]"<<std::endl;
      exit(-1);
   }

	std::ifstream fin(argv[1]);
	if( ! fin ) {
		std::cerr<<"Error(main): cannot read from sac list "<<argv[1]<<std::endl;
		exit(-2);
	}
	bool isfirst=true;
	SacRec sacsum;
	int nadd=0;
	for(std::string line; std::getline(fin, line); nadd++) {
		std::size_t send = line.find_first_of(" ");
		auto line_end = send==std::string::npos ? line.end() : line.begin()+send;
		std::string sacname( line.begin(), line_end );
		if( isfirst ) {
			sacsum.Load(sacname);
			isfirst = false;
		} else {
			SacRec saccur(sacname);
			saccur.Load();
			sacsum.Addf(saccur);
		}
	}
	sacsum.Mul(1./nadd);
	
   sacsum.Write(argv[2]);
   return 0;
}
