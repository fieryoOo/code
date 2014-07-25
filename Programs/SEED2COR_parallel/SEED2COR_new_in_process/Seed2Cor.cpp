#include "InfoLists.h"
#include "SeedRec.h"
#include "SacRec.h"
#include "CCDatabase.h"
#include <iostream>

int main(int argc, char *argv[]) {

   if(argc!=2) {
      std::cerr<<"Usage: "<<argv[0]<<" [Parameters_file]"<<std::endl;
      return 0;
   }

   /* Initialize the CC Database with the input parameter file */
   CCDatabase cdb( argv[1] );

   const CCPARAM cdbParams = cdb.GetParams();
   
	/* iterate through the database and handle all possible events */
	for(DailyInfo dinfo; cdb.GetRec(dinfo); cdb.NextRec()) {
		/* daily info from the database */
		std::cerr<<dinfo.seedname<<" "<<dinfo.staname<<" "<<dinfo.chname<<std::endl;

		/* extract the original sac from seed */
		float gapfrac;
		SacRec sac;
		SeedRec seedcur( dinfo.seedname.c_str(), dinfo.rdsexe.c_str() );
		if( ! seedcur.ExtractSac( dinfo.staname, dinfo.chname, dinfo.sps, dinfo.rec_outname,
										  dinfo.resp_outname, gapfrac, sac ) )	continue;
		sac.Write( dinfo.osac_outname.c_str() );

		/* remove response and cut */
		sac.RmRESP( dinfo.resp_outname.c_str(), dinfo.perl*0.76923, dinfo.perh*1.42857 );
		char evtime[15];
		sprintf( evtime, "%04d%02d%02d000000\0", dinfo.year, dinfo.month, dinfo.day );
		sac.ZoomToEvent( evtime, -12345., -12345., dinfo.t1, dinfo.tlen );
		sac.Write( dinfo.fsac_outname.c_str() );
sac.WriteHD("/usr/temp.SAC");

		/* apply normalizations and convert to am/ph */
		
	}

   return 0;
}
