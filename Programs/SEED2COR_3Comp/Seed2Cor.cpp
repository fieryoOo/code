#define DEBUG
#define MAIN
#include "InfoLists.h"
#include "SeedRec.h"
#include "SacRec.h"
#include "CCDatabase.h"
#include "MyLogger.h"
#include "MyOMP.h"
#include "SysTools.h"
#include "CCRec.h"
#include "StaPair.h"
#include "Timer.h"
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unistd.h>

/* normalize all sac files in sacV (simultaneously if SyncNorm==true) by Ye Tian*/
void TNormAll( std::deque<SacRec>& sacV, const std::vector<DailyInfo>& dinfoV, bool SyncNorm);
/* normalize and apply taper by Ye Tian*/
void FNormAll( std::deque<SacRec>& sacV, const std::vector<DailyInfo>& dinfoV, bool SyncNorm);
bool FileExists(const char* filename);
bool FileExists(const std::string & filename);
std::vector < std::string > GchannelList(CCDatabase & CDB);
/* Generate CC to do List, by Lili Feng */
void GCCtodoList( std::vector<StaPair> & CC_List, std::vector < StaInfo > & StaList, std::vector < std::string > & MDIR, std::vector < std::vector < std::string > > & CHAll, int skipccflag);
/* convert amp & ph files to CC, by Lili Feng */
void CCList2CC( std::vector<StaPair> & CC_List, std::vector < std::vector < std::string > > ChanAll, CCDatabase & CDB);
/* daily component correction, by Lili Feng*/
void Rotation_daily(std::deque < SacRec > & SACV, std::vector < std::string > Rout);
extern MyLogger logger;
MyLogger logger;
extern MEMO memo;
MEMO memo;

class RESPList {
public:
	void Load(const std::vector<std::string> &fileV) {
		_hashT.clear();
		for( auto &fname : fileV ) {
			std::string ch = fname.substr(fname.find_last_of('.')+1);
			_hashT.emplace(ch, fname);
		}
	}
	bool Get(const std::string& ch, std::string& str) {
		try { str = _hashT.at(ch); return true; } 
		catch(const std::out_of_range& oor) { return false; }
	}
private:
	std::unordered_map<std::string, std::string> _hashT;
};

int main(int argc, char *argv[]) {
	if(argc!=2) {
		std::cerr<<"Usage: "<<argv[0]<<" [Parameters_file]"<<std::endl;
		return -1;
	}
	bool SyncNorm = true;
	logger.Rename("logs_Seed2Cor/log_Seed2Cor");

	try {
		/* Initialize the CC Database with the input parameter file */
		CCDatabase cdb( argv[1] );
		const CCPARAM& cdbP = cdb.GetParams();
		/* check total memory available */
		float MemTotal = memo.MemTotal();
		logger.Hold( INFO, "Estimated total memory = "+std::to_string(MemTotal)+" Mb", FuncName );
		logger.flush();
		/* iterate through the database and handle all possible events */
		#pragma omp parallel
		{ // parallel region S
		while( true ) {	// main loop
Timer timer;
std::stringstream sst;
			//int ithread = omp_get_thread_num();
			/* dynamically assign events to threads, one at a time */
			std::vector<DailyInfo> dinfoV;
			bool got;
			#pragma omp critical(cdb)
			{// critical S
			got = cdb.GetRec_AllCH(dinfoV);	// all channels for the current event and station
			cdb.NextEvent();
			} // critical E
			if( !got ) break;

			if( cdbP.fskipesac==2 && cdbP.fskipresp==2 && cdbP.fskipamph==2 ) break;

sst<<dinfoV[0].staname<<"(t0)="<<timer.SecondsElapsed()<<" ";
RESPList respList;
if( cdbP.fskipesac==1 && !dinfoV.empty() ) {
	std::vector<std::string> resp_list;
	List( dinfoV[0].outdir.c_str(), ("RESP."+dinfoV[0].ntname+"."+dinfoV[0].staname+".*").c_str(), 2, resp_list);
	respList.Load(resp_list);
sst<<dinfoV[0].staname<<"(t1)="<<timer.SecondsElapsed()<<" ";
}

			try {	// handle current event
				std::deque<SacRec> sacV;
				std::vector < std::string > sacfout_R;
				std::vector < std::string > sacfout_O;
				std::vector<std::stringstream> reportV( dinfoV.size() );
				std::vector < std::string > ExtractedStaLst, RemovedStaLst, ExistingStaLst;
				// loop through channel list
				for( int ich=0; ich<dinfoV.size(); ich++ ) {
					auto& dinfo = dinfoV[ich];
sst<<dinfo.chname<<"(t2)="<<timer.SecondsElapsed()<<" ";
					bool extract_flag=false;
					/* daily info from the database */
					//logger.Hold( INFO, dinfo.seedname + " " + dinfo.staname + " " + dinfo.chname, FuncName );
					/* stringstream for reporting */
					auto& report = reportV[ich];

					/* extract the original sac from seed */
					SacRec sac( report );
					if ( cdbP.fskipesac==0 || ( ! FileExists(dinfo.osac_outname) && cdbP.fskipesac==1 ) ) {
						/// 1.02-->1.03, check if SAC files are extracted
						float gapfrac;
						sac.SetMaxMemForParallel( MemTotal * dinfo.memomax * 0.8 / omp_get_num_threads() );
						SeedRec seedcur( dinfo.seedname, dinfo.rdsexe, report );
						if( seedcur.ExtractSac( dinfo.staname, dinfo.ntname, dinfo.chname, dinfo.sps, dinfo.rec_outname,
														dinfo.resp_outname, gapfrac, sac ) ) {
							extract_flag=true;
							sac.Write( dinfo.osac_outname );
							ExtractedStaLst.push_back(dinfo.staname+"_"+dinfo.chname);
						}
						sacfout_O.push_back(dinfo.osac_outname);
					} 

					// remove response and cut
					if (  cdbP.fskipresp==2 || (FileExists(dinfo.fsac_outname) && cdbP.fskipresp==1)  ) {
						/// 1.02-->1.03, check if removed response SAC file exists or not
						if (FileExists(dinfo.fsac_outname)) {
							sac.Load(dinfo.fsac_outname);
							sacfout_R.push_back(dinfo.fsac_outname);
						}
						sacfout_O.push_back(dinfo.osac_outname);
						sacV.push_back( std::move(sac) );
						continue;
					}
					if (! extract_flag) {	/// !!! TODO!!! NEED to ensure we get the right resp file!
/*
						std::vector<std::string> resp_list;
						if (List((dinfo.outdir).c_str(), ("RESP."+dinfo.ntname+"."+dinfo.staname+".*."+dinfo.chname).c_str(), 2, resp_list)) {
							dinfo.resp_outname=resp_list[0];
						} else {
*/
						if( ! respList.Get(dinfo.chname, dinfo.resp_outname) ) {
							if (FileExists(dinfo.fsac_outname))	{
								sac.Load(dinfo.fsac_outname);
								sacfout_R.push_back(dinfo.fsac_outname);
							}
							sacfout_O.push_back(dinfo.osac_outname);
							sacV.push_back( std::move(sac) );
							continue;
						}
						if (FileExists(dinfo.osac_outname))	{
							//logger.Hold( INFO,"Reading exsting SAC file: "+dinfo.osac_outname, FuncName );
							//std::cout<<"Reading exsting SAC file: "<<dinfo.osac_outname<<std::endl;
							ExistingStaLst.push_back(dinfo.staname+"_"+dinfo.chname);
							sac.Load(dinfo.osac_outname);
							sacfout_O.push_back(dinfo.osac_outname);
						} else {
							if (FileExists(dinfo.fsac_outname))	{
								sac.Load(dinfo.fsac_outname);
								sacfout_R.push_back(dinfo.fsac_outname);
							}
							sacfout_O.push_back(dinfo.osac_outname);
							sacV.push_back( std::move(sac) );
							continue;
						}
					}
					if( std::min(sac.shd.e, dinfo.tlen+dinfo.t1-sac.shd.b) - std::max(sac.shd.b, dinfo.t1) < cdbP.mintlen ) continue;
					//if ( sac.shd.npts <=1 ) continue;	// ?????

					// remove response, output velocity
					sac.RmRESP( dinfo.resp_outname, dinfo.perl*0.8, dinfo.perh*1.3, dinfo.evrexe, 1 );
					char evtime[15];
					sprintf( evtime, "%04d%02d%02d000000\0", dinfo.year, dinfo.month, dinfo.day );
					sac.ZoomToEvent( evtime, -12345., -12345., dinfo.t1, dinfo.tlen );
					sac.Write( dinfo.fsac_outname );
					// logger.Hold( INFO,"Removed resp SAC file: "+dinfo.fsac_outname, FuncName );
					RemovedStaLst.push_back(dinfo.staname+"_"+dinfo.chname);
					//std::cout<<"Removed resp SAC file: "<<dinfo.fsac_outname<<std::endl;
					//sac.WriteHD("/usr/temp.SAC");
					sacV.push_back( std::move(sac) );
					sacfout_R.push_back(dinfo.fsac_outname);
sst<<dinfo.chname<<"(t3)="<<timer.SecondsElapsed()<<" ";
				}

sst<<dinfoV[0].chname<<"(t4)="<<timer.SecondsElapsed()<<" ";
				// log results for the day
				std::stringstream ssdate;
				ssdate << dinfoV[0].year<<"."<<dinfoV[0].month<<"."<<dinfoV[0].day;
				std::string date(ssdate.str());
				logger.Hold( INFO, "*** Pre-processing done for station "+dinfoV[0].staname+" on "+date+" ***");
				auto LogAllSta = [&]( std::string msg, std::vector<std::string>& sV ) {
					if( sV.empty() ) return;
					msg += " :";
					for(const auto& sta : sV ) msg += " " + sta;
					logger.Hold( INFO, msg, FuncName ); sV.clear();
				};
sst<<dinfoV[0].chname<<"(t5)="<<timer.SecondsElapsed()<<" ";
				LogAllSta( "SAC extraction", ExtractedStaLst );
				LogAllSta( "File existed", ExistingStaLst );
				LogAllSta( "RESP removal", RemovedStaLst );
sst<<dinfoV[0].chname<<"(t6)="<<timer.SecondsElapsed()<<" ";
logger.Hold( INFO, " debug timer: "+sst.str() );
				logger.flush();
logger.Hold( INFO, " debug timer after flush! "+std::to_string(timer.SecondsElapsed()) );

if ( cdbP.fskipamph==2 ) continue;
				if ( cdbP.fskipamph==2 || ( cdbP.fskipamph==1 && 
						FileExists(dinfoV[0].fsac_outname+".am")  && FileExists(dinfoV[0].fsac_outname+".ph")  &&
						FileExists(dinfoV[1].fsac_outname+".am")  && FileExists(dinfoV[1].fsac_outname+".ph") &&
						FileExists(dinfoV[2].fsac_outname+".am")  && FileExists(dinfoV[2].fsac_outname+".ph")  ) ) continue;
std::cerr<<" Error!!! passed!"<<std::endl;

				/* time-domain normalization */
				TNormAll( sacV, dinfoV, SyncNorm );

				/* fre-domain normalization */
				// convert sacs to am&ph and store amp in sacV
				for( int isac=0; isac<sacV.size(); isac++ ) {
					auto& dinfo = dinfoV[isac];
					auto& sac = sacV[isac];
					if( ! sac.sig ) continue;
					SacRec sac_am, sac_ph;
					sac.ToAmPh( sac_am, sac_ph );
					sac = std::move(sac_am); ///???
					//sac = sac_am;
					sac_ph.Write( dinfo.fsac_outname + ".ph" );
				}
				// normalize
				FNormAll( sacV, dinfoV, SyncNorm );
				//FNormAll( sacV, dinfoV, false );
				// write am
				for( int isac=0; isac<sacV.size(); isac++ ) {
					auto& dinfo = dinfoV[isac];
					auto& sac = sacV[isac];
					if( ! sac.sig ) continue;
					sac.Write( dinfo.fsac_outname + ".am" );
				}
				// End of removed SAC file to amp&ph file
				if (cdbP.fdelosac==1 || cdbP.fdelosac==3)
					for (int isac=0; isac<sacfout_O.size(); isac++ )
						if (FileExists(sacfout_O[isac])) fRemove(sacfout_O[isac].c_str());/// 1.02-->1.03
				if (cdbP.fdelosac==2 || cdbP.fdelosac==3)
					for (int isac=0; isac<sacfout_R.size(); isac++ )
						if (FileExists(sacfout_R[isac])) fRemove(sacfout_R[isac].c_str());/// 1.02-->1.03
				// log (if there's any) warnings
				for( const auto& report : reportV )	{
					std::string warning = report.str();
					if( ! warning.empty() )	logger.Hold( WARNING, "\n" + warning, FuncName );
				}
				logger.flush();
			} catch ( std::exception& e ) {
				logger.Hold( ERROR, e.what(), FuncName );
			} // current event done
		} // main while loop
		} // parallel region E

		/*---- CC code start here, by Lili Feng----*/
		if(cdbP.fskipcrco == 3) return 0;
		int fskipcc=cdbP.fskipcrco;
		// Get channel list
		std::vector < std::string > chanL1 = cdb.GchannelList( 1 );
		std::vector < std::string > chanL2 = cdb.GchannelList( 2 );
		std::vector < std::string > chanL3 = cdb.GchannelList( 3 );
		std::vector < std::vector < std::string > > ChannelAll;
		if ( ! chanL1.empty() )
			ChannelAll.push_back(chanL1);
		if ( ! chanL2.empty() )
			ChannelAll.push_back(chanL2);
		if ( ! chanL3.empty() )
			ChannelAll.push_back(chanL3);
		std::vector < StaInfo > stationlist=cdb.GStaList();
		std::vector < std::string > monthdir=cdb.GMonLst();
		std::vector< StaPair > CC_todolist;
		GCCtodoList( CC_todolist, stationlist, monthdir,  ChannelAll, fskipcc ); // generate CC_todolist
		CCList2CC( CC_todolist, ChannelAll, cdb ); // Do CC according to CC_todolist
		/*---- CC code end here ----*/
		logger.Hold( INFO, "All threads finished.", FuncName );
	} catch ( std::exception& e ) {
		logger.Hold(FATAL, e.what(), FuncName);
		return -2;
	} catch (...) {
		logger.Hold(FATAL, "unknown exception", FuncName);
		return -2;
	}
	return 0;
}
