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
	for(std::string line; std::getline(fin, line); ) {
		std::string sacname( line.begin(), line.begin()+line.find_first_of(" \n") );
		if( isfirst ) {
			sacsum.Load(sacname);
		} else {
			SacRec saccur(sacname);
			saccur.Load();
			sacsum.Addf(saccur);
		}
	}

   sacsum.Write(argv[2]);
   return 0;
}
