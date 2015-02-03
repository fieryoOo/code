#include <SacRec.h>
#include <iostream>
#include <stdexcept>

int main( int argc, char* argv[] ) {
   if( argc != 4 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [sac_name] [header_field_name] [value]"<<std::endl;
      exit(-1);
   }
	try {
		SacRec sacrec(argv[1]);
		sacrec.LoadHD();
		sacrec.ChHdr(argv[2], argv[3]);
		sacrec.WriteHD(argv[1]);
	} catch ( const std::exception& e ) {
		std::cerr<<"Error(main): change header failed with exception "<<e.what()<<std::endl;
		return -2;
	}
	return 0;
}
