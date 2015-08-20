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
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <unistd.h>

/* normalize all sac files in sacV (simultaneously if SyncNorm==true) by Ye Tian*/
void TNormAll( std::deque<SacRec>& sacV, const std::vector<DailyInfo>& dinfoV, bool SyncNorm);
/* normalize and apply taper by Ye Tian*/
void FNormAll( std::deque<SacRec>& sacV, const std::vector<DailyInfo>& dinfoV, bool SyncNorm);
bool FileExists(const char* filename);
bool FileExists(const std::string & filename);
std::vector < std::string > GchannelList(CCDatabase & CDB);
void GStaMonList( std::vector < Station > & StaList, std::vector < std::string > & MDIR, CCDatabase & CDB);
void GCCtodoList(std::vector < Station > & StaList, std::vector < std::string > & MDIR, std::vector<CC_todo> & CC_List, std::vector < std::string > & CHAN, int skipccflag);
/* convert amp & ph files to CC, by Lili Feng */
void CCList2CC( std::vector<CC_todo> & CC_List, std::vector < std::string > & CHAN, CCDatabase & CDB);
/* daily component correction, by Lili Feng*/
void Rotation_daily(std::deque < SacRec > & SACV, std::vector < std::string > Rout);
/* stack CC data according to CCstack, by Lili Feng*/
void StackAll(std::vector < Station > & StaList, std::vector < std::string > & MDIR, std::vector < std::string > & CHAN, CCDatabase & CDB);
//bool checkflags(int esacflag, int respflag, int ampflag);
/* rotate CC of ENZ to RTZ, by Lili Feng*/
bool Rotation(sta_pairs & stapair, std::vector < std::string > & CHAN, std::vector < std::string > & DIR, std::string CHPRE);
extern MyLogger logger;
MyLogger logger;
extern MEMO memo;
MEMO memo;

int main(int argc, char *argv[])
{
    if(argc!=2)
    {
        std::cerr<<"Usage: "<<argv[0]<<" [Parameters_file]"<<std::endl;
        return -1;
    }
    bool SyncNorm = true;
    logger.Rename("logs_Seed2Cor/Test");

    try
    {
        /* Initialize the CC Database with the input parameter file */
        CCDatabase cdb( argv[1] );
        //const CCPARAM& cdbParams = cdb.GetParams();
        /* check total memory available */
        float MemTotal = memo.MemTotal();
        logger.Hold( INFO, "Estimated total memory = "+std::to_string(MemTotal)+" Mb", FuncName );
        logger.flush();

        /* iterate through the database and handle all possible events */
        #pragma omp parallel
        {
            // parallel region S
            while( 1 )   // main loop
            {
                //int ithread = omp_get_thread_num();
                /* dynamically assign events to threads, one at a time */
                bool got;
                std::vector<DailyInfo> dinfoV;
                #pragma omp critical(cdb)
                {
                    // critical S
                    got = cdb.GetRec_AllCH(dinfoV);
                    cdb.NextEvent();
                } // critical E
                if( !got ) break;
                if( cdb.GetParams().fskipesac==2 && cdb.GetParams().fskipresp==2 && cdb.GetParams().fskipamph==2 )
                    break;
                try   // handle current event
                {
                    std::deque<SacRec> sacV;
                    std::vector < std::string > sacfout_R;
                    std::vector < std::string > sacfout_O;
                    std::vector<std::stringstream> reportV( dinfoV.size() );
                    /*----  STEP 1: seed to SAC file(removed response) ----*/
                    for( int ich=0; ich<dinfoV.size(); ich++ )
                    {
                        auto& dinfo = dinfoV[ich];
                        bool extract_flag=false;
                        /* daily info from the database */
                        logger.Hold( INFO, dinfo.seedname + " " + dinfo.staname + " " + dinfo.chname, FuncName );
                        //std::cout<<"Scanning:"+dinfo.seedname + " " + dinfo.staname + " " + dinfo.chname<<std::endl;
                        /* stringstream for reporting */
                        auto& report = reportV[ich];
                        /* extract the original sac from seed */
                        SacRec sac( report );

                        if ( cdb.GetParams().fskipesac==0 || ( ! FileExists(dinfo.osac_outname) && cdb.GetParams().fskipesac==1 ) ) ///
                        {
                            /// 1.02-->1.03, check if extract SAC file or not
                            float gapfrac;
                            sac.SetMaxMemForParallel( MemTotal * dinfo.memomax * 0.8 / omp_get_num_threads() );
                            SeedRec seedcur( dinfo.seedname, dinfo.rdsexe, report );
                            if( seedcur.ExtractSac( dinfo.staname, dinfo.chname, dinfo.sps, dinfo.rec_outname,
                                                    dinfo.resp_outname, gapfrac, sac ) )
                            {
                                extract_flag=true;
                                sac.Write( dinfo.osac_outname );
                                logger.Hold( INFO,"Extracting SAC file: "+dinfo.osac_outname, FuncName );
                                //std::cout<<"Extracting SAC file: "<<dinfo.osac_outname<<std::endl;
                                sacfout_O.push_back(dinfo.osac_outname);
                            }
                        }

                        /* remove response and cut */
                        if (  cdb.GetParams().fskipresp==2 || (FileExists(dinfo.fsac_outname) && cdb.GetParams().fskipresp==1)  ) // if do not do resp_sac step
                        {
                            /// 1.02-->1.03, check if removed response SAC file exists or not
                            if (FileExists(dinfo.fsac_outname))
                            {
                                sac.Load(dinfo.fsac_outname);
                                sacfout_R.push_back(dinfo.fsac_outname);
                            }
                            sacfout_O.push_back(dinfo.osac_outname);
                            sacV.push_back( std::move(sac) );
                            continue;
                        }
                        //FIND resp file!
                        if (extract_flag==false)  /// !!! TODO!!! NEED to ensure we get the right resp file!
                        {
                            std::vector<std::string> resp_list;
                            if (List((dinfo.outdir).c_str(), ("RESP*."+dinfo.staname+".*."+dinfo.chname).c_str(), 2, resp_list))
                                dinfo.resp_outname=resp_list[0];
                            else
                            {
                                if (FileExists(dinfo.fsac_outname))
                                {
                                    sac.Load(dinfo.fsac_outname);
                                    sacfout_R.push_back(dinfo.fsac_outname);
                                }
                                sacfout_O.push_back(dinfo.osac_outname);
                                sacV.push_back( std::move(sac) );
                                continue;
                            }
                            if (FileExists(dinfo.osac_outname))
                            {
                                logger.Hold( INFO,"Reading exsting SAC file: "+dinfo.osac_outname, FuncName );
                                //std::cout<<"Reading exsting SAC file: "<<dinfo.osac_outname<<std::endl;
                                sac.Load(dinfo.osac_outname);
                                sacfout_O.push_back(dinfo.osac_outname);
                            }
                            else
                            {
                                if (FileExists(dinfo.fsac_outname))
                                {
                                    sac.Load(dinfo.fsac_outname);
                                    sacfout_R.push_back(dinfo.fsac_outname);
                                }
                                sacfout_O.push_back(dinfo.osac_outname);
                                sacV.push_back( std::move(sac) );
                                continue;
                            }
                        }
                        if ( sac.shd.npts <=1 ) { // Added 1.05--->1.06
                                std::cout<<"ATTENTION: Skip Preprocessing for: "<<dinfo.fsac_outname<<std::endl;
                                continue;
                        }

                        sac.RmRESP( dinfo.resp_outname, dinfo.perl*0.8, dinfo.perh*1.3, dinfo.evrexe );
                        //std::cout<<"Extracting SAC file: "<<dinfo.osac_outname<<std::endl;
                        char evtime[15];
                        sprintf( evtime, "%04d%02d%02d000000\0", dinfo.year, dinfo.month, dinfo.day );
                        sac.ZoomToEvent( evtime, -12345., -12345., dinfo.t1, dinfo.tlen );
                        sac.Write( dinfo.fsac_outname );
                        logger.Hold( INFO,"Removed resp SAC file: "+dinfo.fsac_outname, FuncName );
                        //std::cout<<"Removed resp SAC file: "<<dinfo.fsac_outname<<std::endl;
                        //sac.WriteHD("/usr/temp.SAC");
                        sacV.push_back( std::move(sac) );
                        sacfout_R.push_back(dinfo.fsac_outname);
                    }
                    /*--- Component correction ---*/
                    /*float del=3.0;
                    if (sacV.size()==dinfoV.size() && sacV.size()> 1 &&sacV.size()<4)
                    {
                        // Do not take do esac, not do rsac, do amp into consideration!!!
                        if ( abs(sacV[0].shd.cmpaz-90.0 ) > del && sacV[0].sig && sacV[1].sig)
                            Rotation_daily(sacV, sacfout_R);
                    } */
                    if ( cdb.GetParams().fskipamph==2 ||
                            (  FileExists(dinfoV[0].fsac_outname+".am")  && FileExists(dinfoV[0].fsac_outname+".ph")  &&
                               FileExists(dinfoV[1].fsac_outname+".am")  && FileExists(dinfoV[1].fsac_outname+".ph") &&
                               FileExists(dinfoV[2].fsac_outname+".am")  && FileExists(dinfoV[2].fsac_outname+".ph")  )
                            && cdb.GetParams().fskipamph ==1) continue; /// 1.02-->1.03, check if convert SAC to amp&ph file

                    /* time-domain normalization */
                    TNormAll( sacV, dinfoV, SyncNorm );

                    /* fre-domain normalization */
                    // convert sacs to am&ph and store amp in sacV
                    for( int isac=0; isac<sacV.size(); isac++ )
                    {
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
                    // write am
                    for( int isac=0; isac<sacV.size(); isac++ )
                    {
                        auto& dinfo = dinfoV[isac];
                        auto& sac = sacV[isac];
                        if( ! sac.sig ) continue;
                        sac.Write( dinfo.fsac_outname + ".am" );
                    }
                    // End of removed SAC file to amp&ph file
                    if (cdb.GetParams().fdelosac==1)
                        for (int isac=0; isac<sacfout_O.size(); isac++ )
                            if (FileExists(sacfout_O[isac])) fRemove(sacfout_O[isac].c_str());/// 1.02-->1.03
                    if (cdb.GetParams().fdelosac==2)
                        for (int isac=0; isac<sacfout_R.size(); isac++ )
                            if (FileExists(sacfout_R[isac])) fRemove(sacfout_R[isac].c_str());/// 1.02-->1.03
                    if (cdb.GetParams().fdelosac==3)
                    {
                        for (int isac=0; isac<sacfout_R.size(); isac++ )
                            if (FileExists(sacfout_R[isac])) fRemove(sacfout_R[isac].c_str()); /// 1.02-->1.03
                        for (int isac=0; isac<sacfout_O.size(); isac++ )
                            if (FileExists(sacfout_O[isac])) fRemove(sacfout_O[isac].c_str());/// 1.02-->1.03
                    }
                    /* log if any warning */
                    for( const auto& report : reportV )
                    {
                        std::string warning = report.str();
                        if( ! warning.empty() )
                            logger.Hold( WARNING, "\n" + warning, FuncName );
                        logger.flush();
                    } // for dinfo
                }
                catch ( std::exception& e )
                {
                    logger.Hold( ERROR, e.what(), FuncName );
                } // current event done
            } // main while loop
        } // parallel region E

        /*---- CC code start here, by Lili Feng----*/
        clock_t time_before;
        time_before = clock();
        if(cdb.GetParams().fskipcrco == 3) return 0;
        int fskipcc=cdb.GetParams().fskipcrco;
        std::vector < std::string > channel = GchannelList( cdb ); // generate channel list
        std::string CH_pre;
        CH_pre.append(channel[0].begin(),channel[0].end()-1);
        std::vector < Station > stationlist;
        std::vector < std::string > monthdir;
        std::vector< CC_todo > CC_todolist;
        GStaMonList( stationlist, monthdir, cdb ); // generate station list and monthdir list
        GCCtodoList( stationlist, monthdir, CC_todolist, channel, fskipcc ); // generate CC_todolist
        CCList2CC( CC_todolist, channel, cdb ); // Do CC according to CC_todolist
        StackAll( stationlist, monthdir, channel, cdb );
        /*---- Do Rotation ----*/
        std::vector <sta_pairs> stapair_list;
        for (int s1=0; s1<stationlist.size(); s1++)  // Generate CCtodo List
            for (int s2=0; s2<stationlist.size(); s2++)
            {
                if (stationlist[s1].sta>stationlist[s2].sta || !stationlist[s1].checkdoCC(stationlist[s2]) ) continue;
                stapair_list.push_back(sta_pairs(stationlist[s1].sta,stationlist[s2].sta));
            }
        std::vector< std::string > all_dir;
        all_dir.push_back("COR");
        #pragma omp parallel for
        for (int i=0; i<stapair_list.size(); i++)
            Rotation(stapair_list[i], channel, all_dir, CH_pre);
        /*---- CC code end here ----*/
        logger.Hold( INFO, "All threads finished.", FuncName );
        std::cout<<"Elapsed Time: "<<(float(clock()-time_before))/CLOCKS_PER_SEC<<" secs"<<std::endl;
    }
    catch ( std::exception& e )
    {
        logger.Hold(FATAL, e.what(), FuncName);
        return -2;
    }
    catch (...)
    {
        logger.Hold(FATAL, "unknown exception", FuncName);
        return -2;
    }
    return 0;
}
