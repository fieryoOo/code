#define DEBUG 
#define MAIN

#include "InfoLists.h"
#include "SeedRec.h"
#include "SacRec.h"
#include "CCDatabase.h"
#include "MyLogger.h"
#include "MyOMP.h"
#include <iostream>
#include <sstream>


extern MyLogger logger;
MyLogger logger;

int main(int argc, char *argv[]) {

   if(argc!=2) {
      std::cerr<<"Usage: "<<argv[0]<<" [Parameters_file]"<<std::endl;
      return -1;
   }

	logger.Rename("logs_Seed2Cor/Test");

	try {

		/* Initialize the CC Database with the input parameter file */
		CCDatabase cdb( argv[1] );

		//const CCPARAM& cdbParams = cdb.GetParams();

		/* iterate through the database and handle all possible events */
		#pragma omp parallel
		{ // parallel region S
		while( 1 ) {
			/* dynamically assign events to threads, one at a time */
			bool got;
			DailyInfo dinfo;
			#pragma omp critical
			{ // critical S
			got = cdb.GetRec(dinfo);
			cdb.NextRec();
			} // critical E
			if( !got ) break;

			/* daily info from the database */
			logger.Hold( INFO, dinfo.seedname + " " + dinfo.staname + " " + dinfo.chname, FuncName );

			try {
				/* stringstream for reporting */
				std::stringstream report;

				/* extract the original sac from seed */
				float gapfrac;
				SacRec sac( report );
				SeedRec seedcur( dinfo.seedname, dinfo.rdsexe, report );
				if( ! seedcur.ExtractSac( dinfo.staname, dinfo.chname, dinfo.sps, dinfo.rec_outname,
							dinfo.resp_outname, gapfrac, sac ) )	continue;
				sac.Write( dinfo.osac_outname );

				/* remove response and cut */
				sac.RmRESP( dinfo.resp_outname, dinfo.perl*0.76923, dinfo.perh*1.42857 );
				char evtime[15];
				sprintf( evtime, "%04d%02d%02d000000\0", dinfo.year, dinfo.month, dinfo.day );
				sac.ZoomToEvent( evtime, -12345., -12345., dinfo.t1, dinfo.tlen );
				sac.Write( dinfo.fsac_outname );
				//sac.WriteHD("/usr/temp.SAC");

				/* apply normalizations and convert to am/ph */
				switch( dinfo.tnorm_flag ) {
					case 1: sac.OneBit(); 
						break;
					case 2: sac.RunAvg( dinfo.timehlen, dinfo.Eperl, dinfo.Eperh ); 
						break;
					case 3: throw std::runtime_error("EqkCut not exist!");
						break;
				}
				//sac.Write( "Test_RunAvg.sac" );

				/* convert to freq domain */
				SacRec sac_am, sac_ph;
				sac.ToAmPh( sac_am, sac_ph );
				/* whitening */
				sac_am.RunAvg( dinfo.frechlen, -1., -1. );
				float fl=1./dinfo.perh, fh=1./dinfo.perl;
				sac_am.cosTaperL( fl*0.8, fl );
				sac_am.cosTaperR( fh, fh*1.2 );
				sac_am.Write( dinfo.fsac_outname + ".am" );
				sac_ph.Write( dinfo.fsac_outname + ".ph" );

				/* log if any warning */
				std::string warning = report.str();
				if( ! warning.empty() )
					logger.Hold( WARNING, "\n" + warning, FuncName );
			} catch ( std::exception& e ) {
				logger.Hold( ERROR, e.what(), FuncName );
			}

			logger.flush();
		} // while
		} // parallel region E

	} catch ( std::exception& e ) {
		logger.Hold(FATAL, e.what(), FuncName);
		return -2;
	} catch (...) {
		logger.Hold(FATAL, "unknown exception", FuncName);
		return -2;
	}

	return 0;
}
