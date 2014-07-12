#include "SeedStaInfo.h"
#include "SeedRec.h"
#include "SacRec.h"
#include "CCDatabase.h"
#include <iostream>

int main(int argc, char *argv[]) {

   if(argc!=2) {
      std::cerr<<"Usage: "<<argv[0]<<" [Parameters_file]"<<std::endl;
      return 0;
   }

   /* test time */
/*
   for(int i=0; i<0; i++) {
      CCDatabase* cdbtmp = new CCDatabase(argv[1]);
      delete cdbtmp;
   }
*/

   /* Initialize the CC Database with the input parameter file */
   CCDatabase cdb( argv[1] );

   const CCPARAM cdbParams = cdb.GetParams();
   std::cerr<<cdbParams.sps<<std::endl;
   
   //cdb.NextRecTest(); // testing the function
   cdb.GetRec();
   while( true ) { 
		DailyInfo dinfo = cdb.GetRec();
std::cerr<<dinfo.seedname<<std::endl;
		SeedRec seedcur( dinfo.seedname.c_str(), dinfo.rdsexe.c_str() );
		float gapfrac;
		SacRec sac;
		if( seedcur.ExtractSac( dinfo, gapfrac, sac ) ) {
			sac.Write( "TEST/TEST.SAC" );
			sac.RmRESP( dinfo.resp_outname.c_str(), 5., 80. );
			sac.ZoomToEvent( "20120901000000", -12345., -12345., 1000., 83000. );
			sac.Write( "TEST/ft_TEST.SAC" );
		}
		
		if( ! cdb.NextRec() ) break;
		
exit(-1);
	}

   return 0;
}
