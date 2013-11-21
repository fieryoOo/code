#include "CCDatabase.h"
#include <iostream>

int main(int argc, char *argv[]) {

   if(argc!=2) {
      std::cerr<<"Usage: "<<argv[0]<<" [Parameters_file]"<<std::endl;
      return 0;
   }

   /* test time */
   for(int i=0; i<0; i++) {
      CCDatabase* cdbtmp = new CCDatabase(argv[1]);
      delete cdbtmp;
   }

   /* Initialize the CC Database with the input parameter file */
   CCDatabase cdb( argv[1] );

   const CCPARAM cdbParams = cdb.GetParams();
   std::cerr<<cdbParams.sps<<std::endl;
   

   cdb.NextRec();

   return 0;
}
