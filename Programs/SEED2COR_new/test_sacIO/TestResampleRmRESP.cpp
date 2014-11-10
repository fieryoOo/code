/* testing the methods of the SacRec object */

#include "SacRec.h"
#include <string>

int main( int argc, char* argv[] )
{
   if( argc != 3 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [input_sac] [RESP file]"<<std::endl;
      exit(-1);
   }
   /* read the sac file into sacrec1 */
   SacRec sacrec(argv[1]);
   sacrec.Load();

   /* remove resp */
   float perl = 0.5, perh = 100.;
   std::string evrexe("/home/yeti4009/Softwares/evalresp-3.3.3/evalresp");
   if( ! sacrec.RmRESP(evrexe.c_str(), argv[2], perl, perh) ) {
      std::cerr<<"RmRESP failed!"<<std::endl;
      return -1;
   }

   std::string outname = "ft_" + std::string(argv[1]);
   sacrec.Write(outname.c_str());

   return 0;
}
