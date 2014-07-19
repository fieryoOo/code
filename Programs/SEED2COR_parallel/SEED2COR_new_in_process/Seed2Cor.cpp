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
   
   //cdb.NextRecTest(); // testing the function
   cdb.GetRec();
   while( true ) { 
		DailyInfo dinfo = cdb.GetRec();
		std::cerr<<dinfo.seedname<<" "<<dinfo.staname<<" "<<dinfo.chname<<std::endl;
/*
		dinfo.chname = "LHZ";
		std::string osac_name = "2012.SEP.1." + dinfo.staname + "." + dinfo.chname + ".SAC";
		dinfo.sps = cdbParams.sps;
		dinfo.rec_outname = osac_name + "_rec1";
		dinfo.resp_outname = "RESP.TA." + dinfo.staname + ".." + dinfo.chname;
*/
		SeedRec seedcur( dinfo.seedname.c_str(), dinfo.rdsexe.c_str() );
		float gapfrac;
		SacRec sac;
		//if( seedcur.ExtractSac( dinfo, gapfrac, sac ) ) {
		if( seedcur.ExtractSac( dinfo.staname, dinfo.chname, dinfo.sps, dinfo.rec_outname, dinfo.resp_outname, gapfrac, sac ) ) {
			sac.Write( dinfo.osac_outname.c_str() );
			sac.RmRESP( dinfo.resp_outname.c_str(), dinfo.perl*0.76923, dinfo.perh*1.42857 );
			sac.ZoomToEvent( "20120901000000", -12345., -12345., dinfo.t1, dinfo.tlen );
			sac.Write( dinfo.fsac_outname.c_str() );
		}
		
		if( ! cdb.NextRec() ) break;
		
//exit(-1);
	}

   return 0;
}
