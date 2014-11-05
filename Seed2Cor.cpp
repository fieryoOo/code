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


void TNormAll( std::vector<SacRec>& sacV, const std::vector<DailyInfo>& dinfoV, bool SyncNorm ) {
	if( sacV.empty() ) return;
	if( sacV.size()!=dinfoV.size() )
		throw std::runtime_error("Error(TNormAll): size mismatch ("+std::to_string(sacV.size())+" - "+std::to_string(dinfoV.size()));

	SacRec sac_sigmax;
	for( int isac=0; isac<sacV.size(); isac++ ) {
		auto& dinfo = dinfoV[isac];
		auto& sac = sacV[isac];
		if( !(sac.sig) ) continue;
		/* apply normalizations and convert to am/ph */
		switch( dinfo.tnorm_flag ) {
			case 1: 
				sac.OneBit(); 
				break;
			case 2: 
				if( SyncNorm ) {
					SacRec sac_sm;
					/* filter into the earthquake band */
					if( dinfo.Eperl != -1. ) {
						SacRec sac_eqk;
						float f2 = 1./dinfo.Eperh, f1 = f2*0.6, f3 = 1./dinfo.Eperl, f4 = f3*1.4;
						sac.Filter( f1, f2, f3, f4, sac_eqk );
						sac_eqk.Smooth( dinfo.timehlen, sac_sm );
					} else {
						sac.Smooth( dinfo.timehlen, sac_sm );
					}
					/* max smoothed signal */
					sac_sigmax.PullUpTo( sac_sm );
				} else {
					sac.RunAvg( dinfo.timehlen, dinfo.Eperl, dinfo.Eperh ); 
				}
				break;
			case 3: 
				throw std::runtime_error("EqkCut not exist!");
				break;
		}
	}

	if( SyncNorm )
		for( auto& sac : sacV )
			if( sac.sig ) sac.Divf( sac_sigmax );
	sacV[0].Write( "Test_RunAvg.sac" );
}


extern MyLogger logger;
MyLogger logger;

int main(int argc, char *argv[]) {

   if(argc!=2) {
      std::cerr<<"Usage: "<<argv[0]<<" [Parameters_file]"<<std::endl;
      return -1;
   }

	bool SyncNorm = true;
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
			std::vector<DailyInfo> dinfoV;
			#pragma omp critical
			{ // critical S
			got = cdb.GetRec_AllCH(dinfoV);
			cdb.NextEvent();
			} // critical E
			if( !got ) break;

			try {
				std::vector<SacRec> sacV;
				std::vector<std::stringstream> reportV( dinfoV.size() );
				/* seed to fsac */
				for( int ich=0; ich<dinfoV.size(); ich++ ) {
					auto& dinfo = dinfoV[ich];
					/* daily info from the database */
					logger.Hold( INFO, dinfo.seedname + " " + dinfo.staname + " " + dinfo.chname, FuncName );
					/* stringstream for reporting */
					auto& report = reportV[ich];

					/* extract the original sac from seed */
					float gapfrac;
					SacRec sac( report );
					SeedRec seedcur( dinfo.seedname, dinfo.rdsexe, report );
					if( ! seedcur.ExtractSac( dinfo.staname, dinfo.chname, dinfo.sps, dinfo.rec_outname,
								dinfo.resp_outname, gapfrac, sac ) ) {
						sacV.push_back( std::move(sac) );
						//sacV.push_back( SacRec() );
						continue;
					}
					sac.Write( dinfo.osac_outname );

					/* remove response and cut */
					sac.RmRESP( dinfo.resp_outname, dinfo.perl*0.76923, dinfo.perh*1.42857 );
					char evtime[15];
					sprintf( evtime, "%04d%02d%02d000000\0", dinfo.year, dinfo.month, dinfo.day );
					sac.ZoomToEvent( evtime, -12345., -12345., dinfo.t1, dinfo.tlen );
					sac.Write( dinfo.fsac_outname );
					//sac.WriteHD("/usr/temp.SAC");
					sacV.push_back( std::move(sac) );
				}

				/* time-domain normalization */
				TNormAll( sacV, dinfoV, SyncNorm );
			
				/* fre-domain normalization */
				for( int isac=0; isac<sacV.size(); isac++ ) {
					auto& dinfo = dinfoV[isac];
					auto& sac = sacV[isac];
					if( ! sac.sig ) continue;
					/* convert to freq domain */
					SacRec sac_am, sac_ph;
					sac.ToAmPh( sac_am, sac_ph );
					/* whitening */
					//sac_am.RunAvg( dinfo.frechlen, -1., -1. );
					float fl=1./dinfo.perh, fh=1./dinfo.perl;
					sac_am.cosTaperL( fl*0.8, fl );
					sac_am.cosTaperR( fh, fh*1.2 );
					sac_am.Write( dinfo.fsac_outname + ".am" );
					sac_ph.Write( dinfo.fsac_outname + ".ph" );

					/* log if any warning */
					std::string warning = reportV.at(isac).str();
					if( ! warning.empty() )
						logger.Hold( WARNING, "\n" + warning, FuncName );

					logger.flush();
				} // for dinfo
			} catch ( std::exception& e ) {
				logger.Hold( ERROR, e.what(), FuncName );
			}
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
