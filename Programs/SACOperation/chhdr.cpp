#include <SacRec.h>
#include <iostream>

int main( int argc, char* argv[] ) {
   if( argc != 4 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [sac_name] [header_field_name] [value]"<<std::endl;
      exit(-1);
   }
   SacRec sacrec(argv[1]);
   sacrec.LoadHD();
   if( ! sacrec.ChHdr(argv[2], argv[3]) ) {
      std::cerr<<"Error(main): change header failed!"<<std::endl;
      exit(0);
   }
   sacrec.WriteHD(argv[1]);
   return 0;
}
