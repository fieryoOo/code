#include "CCDatabase.h"
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>

/*-------------------------------------- CCDatabase ----------------------------------------*/
/* Constructor: read in parameters -> fill in seed file list -> fill in station list */
CCDatabase::CCDatabase( const std::string& fname )
{
    CCParams.Load( fname );
    CCParams.CheckAll();
    seedlst.Load( CCParams.seedfname );
    stalst.Load( CCParams.stafname );
    chlst.Load( CCParams.chlst_info );
    FillDInfo();
}

/* fill dinfo with parameters in CCParams */
void CCDatabase::FillDInfo()
{
    dinfo.rdsexe = CCParams.rdsexe;
    dinfo.evrexe = CCParams.evrexe;
    dinfo.sps = CCParams.sps;
    dinfo.perl = CCParams.perl;
    dinfo.perh = CCParams.perh;
    dinfo.t1 = CCParams.t1;
    dinfo.tlen = CCParams.tlen;
    dinfo.tnorm_flag = CCParams.tnorm_flag;
    dinfo.timehlen = CCParams.timehlen;
    dinfo.Eperl = CCParams.Eperl;
    dinfo.Eperh = CCParams.Eperh;
    dinfo.frechlen = CCParams.frechlen;
    dinfo.memomax = CCParams.memomax;
    dinfo_rdy = false;
}

std::vector < std::string > CCDatabase::GchannelList(const int flag)
{
    Channellist CHLst;
    if (flag==1)
        CHLst.Load( CCParams.chlst1_info );
    else if (flag==2)
        CHLst.Load( CCParams.chlst2_info );
    else if (flag==3)
        CHLst.Load( CCParams.chlst3_info );
    else
        CHLst=chlst;
    std::vector < std::string > channel;
    CHLst.Rewind();
    while (1)
    {
        if (CHLst.IsEnded()) break;
        channel.push_back(*CHLst.GetRec());
        CHLst.NextRec();
    }
    return channel;
}

std::vector < StaInfo > CCDatabase::GStaList()
{
    std::vector < StaInfo > stations;
    stalst.Rewind();
    while (1)
    {
        if (stalst.IsEnded()) break;
        stations.push_back(*stalst.GetRec());
        stalst.NextRec();
    }
    return stations;
}
std::vector < std::string > CCDatabase::GMonLst()
{
    seedlst.Rewind();
    std::vector < std::string > MonList;
    int year=-999, month=-999, Nyear, Nmonth;
    bool newMFlag;
    const std::vector<std::string> MName
    {
        "INVALID",
        "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
        "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
    };
    while (1)
    {
        if (seedlst.IsEnded()) break;
        Nyear=(*(seedlst.GetRec())).year;
        Nmonth=(*(seedlst.GetRec())).month;
        if (year!=Nyear || month !=Nmonth)
            newMFlag=true;
        year=Nyear;
        month=Nmonth;
        if (newMFlag==true)
        {
            MonList.push_back(std::to_string(year)+"."+MName[month]);
            newMFlag=false;
        }
        seedlst.NextRec();
    }
    return MonList;
}


/* pull out the next daily record from the database */
bool CCDatabase::NextRecTest()
{
    std::cerr<<*(seedlst.GetRec())<<std::endl;
    if( seedlst.ReLocate( 2012, 1, 26 ) )
    {
        std::cerr<<*(seedlst.GetRec())<<std::endl;
    }
    if( seedlst.NextRec() ) std::cerr<<*(seedlst.GetRec())<<std::endl;
    if( seedlst.ReLocate( 2011, 2, 26 ) )
    {
        std::cerr<<*(seedlst.GetRec())<<std::endl;
    }
    if( seedlst.NotEnded() ) std::cerr<<*(seedlst.GetRec())<<std::endl;
    if( seedlst.ReLocate( 2012, 2, 26 ) )
    {
        std::cerr<<*(seedlst.GetRec())<<std::endl;
    }
    if( seedlst.NotEnded() ) std::cerr<<*(seedlst.GetRec())<<std::endl;

    std::cerr<<*(stalst.GetRec())<<std::endl;
    if( stalst.ReLocate( "J23A" ) )
    {
        std::cerr<<*(stalst.GetRec())<<std::endl;
    }
    if( stalst.ReLocate( "SAO" ) )
    {
        std::cerr<<*(stalst.GetRec())<<std::endl;
    }
    if( stalst.NextRec() ) std::cerr<<*(stalst.GetRec())<<std::endl;
    if( stalst.ReLocate( "CMA" ) )
    {
        std::cerr<<*(stalst.GetRec())<<std::endl;
    }
    if( stalst.NotEnded() ) std::cerr<<*(stalst.GetRec())<<std::endl;

    return true;
}

bool CCDatabase::NextRec()
{
    if( chlst.IsEnded()  ||
            stalst.IsEnded() ||
            seedlst.IsEnded() ) return false;
    dinfo_rdy = false;
    if( chlst.NextRec() ) return true;
    chlst.Rewind();
    if( stalst.NextRec() ) return true;
    stalst.Rewind();
    if( seedlst.NextRec() ) return true;
    return false;
}

bool CCDatabase::NextEvent()
{
    if( stalst.IsEnded() ||
            seedlst.IsEnded() ) return false;
    if( stalst.NextRec() ) return true;
    stalst.Rewind();
    if( seedlst.NextRec() ) return true;
    return false;
}

void CCDatabase::Rewind()
{
    chlst.Rewind();
    stalst.Rewind();
    seedlst.Rewind();
}

bool CCDatabase::GetRec(DailyInfo& dinfoout)
{
    if( chlst.IsEnded()  ||
            stalst.IsEnded() ||
            seedlst.IsEnded() ) return false;
    if( ! dinfo_rdy )
    {
        dinfo.Update( *(seedlst.GetRec()), *(stalst.GetRec()), *(chlst.GetRec()) );
        dinfo_rdy = true;
    }
    dinfoout = dinfo;
    return true;
}

bool CCDatabase::GetRec_AllCH( std::vector<DailyInfo>& dinfoV )
{
    if( stalst.IsEnded() || seedlst.IsEnded() ) return false;
    for( chlst.Rewind(); !chlst.IsEnded(); chlst.NextRec() )
    {
        dinfo.Update( *(seedlst.GetRec()), *(stalst.GetRec()), *(chlst.GetRec()) );
        dinfoV.push_back(dinfo);
    }
    return true;
}


/*-------------------------------------- CCPARAM ----------------------------------------*/
static int isTermi(int deno)
{
    while(deno%2==0) deno *= 0.5;
    while(deno%5==0) deno *= 0.2;
    return deno==1;
}
bool FindInPath( const std::string fname, std::string& absname )
{
    char* pPath = getenv("PATH");
    if( ! pPath ) return false;
    std::stringstream sPath(pPath);
    for(std::string pathcur; std::getline(sPath, pathcur, ':'); )
    {
        std::string testname = pathcur + '/' + fname;
        std::ifstream fin(testname);
        if( fin )
        {
            absname = testname;
            return true;
        }
    }
    return false;
}
/* read in parameters for the CC Database from the inputfile */
void CCPARAM::Load( const std::string fname )
{
    std::ifstream fin(fname);
    if( ! fin ) throw ErrorCD::BadFile(FuncName, fname);
    int nparam = 0;
    for( std::string stmp; std::getline(fin, stmp); )
    {
        if( Set(stmp) == Succeed ) nparam++;
        else nparam++;
    }
    fin.close();

    //(*report)<<"### "<<nparam<<" succed loads from param file "<<fname<<". ###"<<std::endl;
    std::cout<<"### "<<nparam<<" succed loads from param file "<<fname<<". ###"<<std::endl;
}
CCPARAM::SetRet CCPARAM::Set( const std::string& input )
{
    std::stringstream sin(input);
    std::string field;
    if( !(sin >> field) ) return EmptyInfo;
    if( field.at(0) == '#' )
    {
        return Comment;
    }
    else if( field == "rdsexe" )
    {
        if( !(sin >> rdsexe) ) return BadValue;
    }
    else if( field == "evrexe" )
    {
        if( !(sin >> evrexe) ) return BadValue;
    }
    else if( field == "stafname" )
    {
        if( !(sin >> stafname) ) return BadValue;
    }
    else if( field == "seedfname" )
    {
        if( !(sin >> seedfname) ) return BadValue;
    }
    else if( field == "chlst1" )
    {
        if( !std::getline(sin, chlst1_info, '#') ) return BadValue;
        chlst1_info = chlst1_info.substr( chlst1_info.find_first_not_of(" \t") );
    }
    else if( field == "chlst2" )
    {
        if( !std::getline(sin, chlst2_info, '#') ) return BadValue;
        chlst2_info = chlst2_info.substr( chlst2_info.find_first_not_of(" \t") );
    }
    else if( field == "chlst3" )
    {
        if( !std::getline(sin, chlst3_info, '#') ) return BadValue;
        chlst3_info = chlst3_info.substr( chlst3_info.find_first_not_of(" \t") );
    }
    else if( field == "sps" )
    {
        float ftmp;
        if( !(sin >> ftmp) ) return BadValue;
        if( ftmp != static_cast<int>(ftmp) ) return InvalidValue;
        sps = ftmp;
    }
    else if( field == "gapfrac" )
    {
        if( !(sin >> gapfrac) ) return BadValue;
    }
    else if( field == "t1" )
    {
        if( !(sin >> t1) ) return BadValue;
    }
    else if( field == "tlen" )
    {
        if( !(sin >> tlen) ) return BadValue;
    }
    else if( field == "perl" )
    {
        if( !(sin >> perl) ) return BadValue;
    }
    else if( field == "perh" )
    {
        if( !(sin >> perh) ) return BadValue;
    }
    else if( field == "tnorm_flag" )
    {
        if( !(sin >> tnorm_flag) ) return BadValue;
    }
    else if( field == "Eperl" )
    {
        if( !(sin >> Eperl) ) return BadValue;
    }
    else if( field == "Eperh" )
    {
        if( !(sin >> Eperh) ) return BadValue;
    }
    else if( field == "timehlen" )
    {
        if( !(sin >> timehlen) ) return BadValue;
    }
    else if( field == "frechlen")
    {
        if( !(sin >> frechlen) ) return BadValue;
    }
    else if( field == "fwname" )
    {
        if( !(sin >> fwname) ) return BadValue;
    }
    else if( field == "ftlen" )
    {
        if( !(sin >> ftlen) ) return BadValue;
    }
    else if( field == "fprcs" )
    {
        if( !(sin >> fprcs) ) return BadValue;
    }
    else if( field == "memomax")
    {
        if( !(sin >> memomax) ) return BadValue;
    }
    else if( field == "lagtime")
    {
        if( !(sin >> lagtime) ) return BadValue;
    }
    else if( field == "mintlen" )
    {
        if( !(sin >> mintlen) ) return BadValue;
    }
    else if( field == "fdelosac" )
    {
        if( !(sin >> fdelosac) ) return BadValue;
    }
    else if( field == "fdelamph" )
    {
        if( !(sin >> fdelamph) ) return BadValue;
    }
    else if( field == "fskipesac" )
    {
        if( !(sin >> fskipesac) ) return BadValue;
    }
    else if( field == "fskipresp" )
    {
        if( !(sin >> fskipresp) ) return BadValue;
    }
    else if( field == "fskipamph" )
    {
        if( !(sin >> fskipamph) ) return BadValue;
    }
    else if( field == "fskipcrco" )
    {
        if( !(sin >> fskipcrco) ) return BadValue;
    }
    else if( field == "CorOutflag" )
    {
        if( !(sin >> CorOutflag) ) return BadValue;
    }
    else if( field == "fstackall" )
    {
        if( !(sin >> fstackall) ) return BadValue;
    }
    chlst_info=chlst1_info+chlst2_info+chlst3_info;

    return EmptyInfo;
}
bool CCPARAM::CheckAll()
{
    /* check one parameter at a time (29 total) */
    bool ERR = false;
    std::cout<<"---------------------------Checking Input Parameters---------------------------"<<std::endl;
    //1. rdsexe
    std::cout<<"rdseed excutable:\t"<<rdsexe<<std::endl;
    if( rdsexe.empty() || access(rdsexe.c_str(), R_OK)!=0 )
    {
        std::cerr<<"   Warning: invalid path "<<rdsexe<<std::endl;
        if( FindInPath( "rdseed", rdsexe) )
            std::cerr<<"            corrected to "<<rdsexe<<std::endl;
        else
            ERR = true;
    }
    else if( ! rdsexe.find("rdseed") )
    {
        std::cout<<"   Warning: Are you sure this is an rdseed excutable?"<<std::endl;
    }
    //2. evrexe
    std::cout<<"evalresp excutable:\t"<<evrexe<<std::endl;
    if( evrexe.empty() || access(evrexe.c_str(), R_OK)!=0 )
    {
        std::cerr<<"   Warning: invalid path "<<evrexe<<std::endl;
        if( FindInPath( "evalresp", evrexe) )
            std::cerr<<"            corrected to "<<evrexe<<std::endl;
        else
            ERR = true;
    }
    else if( ! evrexe.find("evalresp") )
    {
        std::cout<<"   Warning: Are you sure this is an evalresp excutable?"<<std::endl;
    }
    //3. stafname
    std::cout<<"station list:\t\t"<<stafname<<std::endl;
    if( stafname.empty() )
    {
        std::cerr<<"   Error: empty file name "<<std::endl;
        ERR = true;
    }
    else
    {
        std::ifstream filetmp(stafname);
        if( ! filetmp )
        {
            std::cerr<<"   Error: cannot access file "<<stafname<<std::endl;
            ERR = true;
        }
        else
        {
            std::string buff;
            std::getline(filetmp, buff);
            char stmp[buff.length()];
            float ftmp;
            if( sscanf(buff.c_str(), "%s %f %f", stmp, &ftmp, &ftmp) != 3 )
            {
                std::cerr<<"   Error: incorrect format in file "<<stafname<<"! Shoud be (sta lon lat)"<<std::endl;
                ERR = true;
            }
            filetmp.close();
        }
    }
    //4. seedfname
    std::cout<<"seed list:\t\t"<<seedfname<<std::endl;
    if( seedfname.empty() )
    {
        std::cerr<<"   Error: empty file name "<<std::endl;
        ERR = true;
    }
    else
    {
        std::ifstream filetmp(seedfname);
        if( ! filetmp )
        {
            std::cerr<<"   Error: cannot access file "<<seedfname<<std::endl;
            ERR = true;
        }
        else
        {
            std::string buff;
            std::getline(filetmp, buff);
            char stmp[buff.length()];
            int itmp;
            if( sscanf(buff.c_str(), "%s %d %d %d", stmp, &itmp, &itmp, &itmp) != 4 )
            {
                std::cerr<<"   Error: incorrect format in file "<<seedfname<<"! (path/seedfile yyyy mm dd) is expected"<<std::endl;
                ERR = true;
            }
            filetmp.close();
        }
    }
    //5-1. chlst1_info
    std::cout<<"channel list 1 info:\t"<<chlst1_info<<std::endl;
    if( chlst1_info.empty() )
    {
        std::cerr<<"   Error: empty string "<<std::endl;
        ERR = true;
    }
    //5-2. chlst2_info
    if( ! chlst2_info.empty() )
        std::cout<<"channel list 2 info:\t"<<chlst2_info<<std::endl;
    //5-3. chlst2_info
    if ( ! chlst3_info.empty() )
        std::cout<<"channel list 2 info:\t"<<chlst3_info<<std::endl;

    //6. sps
    std::cout<<"target sampling rate:\t"<<sps<<std::endl;
    if( sps == NaN )
    {
        std::cerr<<"   Error: empty value"<<std::endl;
        ERR = true;
    }
    else if( sps <= 0 )
    {
        std::cerr<<"   Error: a positive integer is expected!"<<std::endl;
        ERR = true;
    }
    else if( !isTermi(sps) )
    {
        std::cerr<<"   Error: 1/"<<sps<<" isn't a terminating decimal, will cause rounding error!"<<std::endl;
        ERR = true;
    }
    //7. gapfrac
    std::cout<<"max gap fraction:\t"<<gapfrac*100<<"%"<<std::endl;
    if( gapfrac == NaN )
    {
        std::cerr<<"   Error: empty value"<<std::endl;
        ERR = true;
    }
    else if( gapfrac<0 || gapfrac>1 )
    {
        std::cerr<<"   Error: a number between 0. - 1. is expected!"<<std::endl;
        ERR = true;
    }
    //8. t1
    std::cout<<"cutting begining:\t"<<t1<<"sec"<<std::endl;
    if( t1 == NaN )
    {
        std::cerr<<"   Error: empty value"<<std::endl;
        ERR = true;
    }
    else if( fabs(t1) > 86400 )
    {
        std::cout<<"   Warning: "<<t1<<"sec exceeds one day."<<std::endl;
    }
    //9. tlen
    std::cout<<"time-rec length:\t"<<tlen<<"sec"<<std::endl;
    if( tlen == NaN )
    {
        std::cerr<<"   Error: empty value"<<std::endl;
        ERR = true;
    }
    else
    {
        float ftmp = t1+tlen;
        if( ftmp>86400 || ftmp<0 )
            std::cout<<"   Warning: ending time '"<<ftmp<<"sec' out of range."<<std::endl;
    }
    //10,11. perl perh
    std::cout<<"signal per band:\t"<<perl<<" - "<<perh<<"sec"<<std::endl;
    if( perl==NaN || perh==NaN )
    {
        std::cerr<<"   Error: empty value"<<std::endl;
        ERR = true;
    }
    else
    {
        if( perl<0 || perh<0 )
        {
            std::cerr<<"   Error: period band can not go below 0."<<std::endl;
            ERR = true;
        }
        if(perl<2./sps)
        {
            std::cerr<<"   Error: "<<perl<<"sec is out of the lower limit at sps "<<sps<<std::endl;
            ERR = true;
        }
        if(perl<5./sps) std::cout<<"   Warning: signal at "<<perl<<"sec might be lost at a sampling rate of "<<sps<<std::endl;
    }
    //12. tnorm_flag
    std::cout<<"t-norm method:\t\t";
    if( tnorm_flag == NaN )
    {
        std::cerr<<"   Error: empty value"<<std::endl;
        ERR = true;
    }
    else if(tnorm_flag==0) std::cout<<"none"<<std::endl;
    else if(tnorm_flag==1) std::cout<<"One-bit"<<std::endl;
    else if(tnorm_flag==2) std::cout<<"Running average"<<std::endl;
    else if(tnorm_flag==3) std::cout<<"Earthquake cutting"<<std::endl;
    else
    {
        std::cerr<<std::endl<<"   Error: Unknow method. integer between 0 - 3 is expected"<<std::endl;
        ERR = true;
    }
    //13,14. Eperl Eperh
    if(tnorm_flag==2 || tnorm_flag==3)
    {
        if( Eperl==NaN || Eperh==NaN )
        {
            std::cerr<<"   Error: empty value"<<std::endl;
            ERR = true;
        }
        else if( Eperl == -1. )
        {
            std::cout<<"Eqk filter:\t\toff"<<std::endl;
        }
        else
        {
            std::cout<<"Eqk filter:\t\t"<<Eperl<<" - "<<Eperh<<"sec"<<std::endl;
            if( Eperl<0 || Eperh<0 )
            {
                std::cerr<<"   Error: period band can not go below 0."<<std::endl;
                ERR = true;
            }
        }
    }
    //15. timehlen
    if(tnorm_flag==2 || tnorm_flag==3)
    {
        std::cout<<"t-len for run-avg:\t"<<timehlen<<std::endl;
        if( timehlen == NaN )
        {
            std::cerr<<"   Error: empty value"<<std::endl;
            ERR = true;
        }
        else if(timehlen<0)
        {
            std::cerr<<"   Error: positive number is expected"<<std::endl;
            ERR = true;
        }
    }
    //16. frechlen
    std::cout<<"t-len for whitening:\t"<<frechlen<<"  ";
    if( frechlen == NaN )
    {
        std::cerr<<"   Error: empty value"<<std::endl;
        ERR = true;
    }
    else if(frechlen==-1)
    {
        std::cout<<"input smoothing file will be used";
    }
    else if(frechlen<0)
    {
        std::cerr<<std::endl<<"   Error: non-negative number is expected";
        ERR = true;
    }
    std::cout<<std::endl;
    //17. fwname
    if(frechlen==-1)
    {
        std::cout<<"spec reshaping file:\t"<<fwname<<std::endl;
        if( fwname.empty() )
        {
            std::cerr<<"   Error: empty file name"<<std::endl;
            ERR = true;
        }
        else if( access(fwname.c_str(), R_OK)!=0 )
        {
            std::cerr<<"   Error: cannot access file "<<fwname<<std::endl;
            ERR = true;
        }
    }
    //18. ftlen
    if( ftlen == NaN )
    {
        std::cerr<<"   Error: empty value"<<std::endl;
        ERR = true;
    }
    else if(ftlen<=0)
    {
        std::cout<<"cor-time-len correction\toff"<<std::endl;
    }
    else
    {
        std::cout<<"cor-time-len correction\ton"<<std::endl; /*if(tnorm_flag==3) ftlen = 2;*/
    }
    //19. fprcs
    if( fprcs == NaN )
    {
        std::cerr<<"   Error: empty value"<<std::endl;
        ERR = true;
    }
    else if(fprcs<=0)
    {
        std::cout<<"prcsr-signal checking\toff"<<std::endl;
    }
    else
    {
        std::cout<<"prcsr-signal checking\ton"<<std::endl;
    }
    //20. memomax
    std::cout<<"memory fraction:\t"<<memomax*100<<"%"<<std::endl;
    if( memomax == NaN )
    {
        std::cerr<<"   Error: empty value"<<std::endl;
        ERR = true;
    }
    else if( memomax<0 || memomax>1 )
    {
        std::cerr<<"   Error: a number between 0. - 1. is expected!"<<std::endl;
        ERR = true;
    }
    //21. lagtime
    std::cout<<"cor lag time:\t\t"<<lagtime<<"sec"<<std::endl;
    if( lagtime == NaN )
    {
        std::cerr<<"   Error: empty value"<<std::endl;
        ERR = true;
    }
    else if(lagtime<0)
    {
        std::cerr<<"   Error: negative lagtime!"<<std::endl;
        ERR = true;
    }
    else if(lagtime>524288)
    {
        std::cout<<"   Warning: lag time exceeds the maximum"<<std::endl;
    }
    //22. mintlen
    std::cout<<"min time len:\t\t"<<mintlen<<"sec"<<std::endl;
    if( mintlen == NaN )
    {
        std::cerr<<"   Error: empty value"<<std::endl;
        ERR = true;
    }
    else if( mintlen>86400 )
    {
        std::cout<<"   Warning: allowed minimum time length larger than a day."<<std::endl;
    }
    //23. fdelosac
    std::cout<<"Delete orig sacs?\t";
    if( fdelosac == 1 )
    {
        std::cout<<"Yes"<<std::endl;
    }
    else
    {
        std::cout<<"No"<<std::endl;
    }
    //24. fdelamph
    std::cout<<"Delete am&ph files?\t";
    if( fdelamph == 1 )
    {
        std::cout<<"Yes"<<std::endl;
    }
    else
    {
        std::cout<<"No"<<std::endl;
    }
    //25. fskipesac
    std::cout<<"ExtractSac()\t\t";
    if(fskipesac==2) std::cout<<"skip"<<std::endl;
    else if(fskipesac==1) std::cout<<"skip if target file exists"<<std::endl;
    else std::cout<<"overwrite if target file exists"<<std::endl;
    //26. fskipresp
    std::cout<<"RmRESP()\t\t";
    if(fskipresp==2) std::cout<<"skip"<<std::endl;
    else if(fskipresp==1) std::cout<<"skip if target file exists"<<std::endl;
    else std::cout<<"overwrite if target file exists"<<std::endl;
    //27. fskipamph
    std::cout<<"TempSpecNorm()\t\t";
    if(fskipamph==2) std::cout<<"skip"<<std::endl;
    else if(fskipamph==1) std::cout<<"skip if target file exists"<<std::endl;
    else std::cout<<"overwrite if target file exists"<<std::endl;
    //28. fskipcrco
    std::cout<<"CrossCorr()\t\t";
    if(fskipcrco==2) std::cout<<"skip"<<std::endl;
    else if(fskipcrco==1) std::cout<<"skip if target file exists"<<std::endl;
    else std::cout<<"overwrite if target file exists"<<std::endl;
    //29. CorOutflag
    std::cout<<"CC Output: \t\t";
    switch( CorOutflag )
    {
    case 0:
        std::cout<<"monthly"<<std::endl;
        break;
    case 1:
        std::cout<<"daily"<<std::endl;
        break;
    case 2:
        std::cout<<"monthly + daily"<<std::endl;
        break;
    default:
        std::cerr<<std::endl<<"   Error: unknow outflag("<<CorOutflag<<")!"<<std::endl;
        ERR = true;
    }
    // 30. StackAll
    if (fstackall==1) std::cout<<"Stack all data."<<std::endl;
    else std::cout<<"Skip stack all"<<std::endl;
    //check for EOF
    std::cout<<"-------------------------------Checking completed-------------------------------"<<std::endl;
    if(ERR) exit(0);
    /* prompt to continue */
    char cin;
    std::string line;
    std::cout<<"Continue?  ";
    std::getline(std::cin, line);
    sscanf(line.c_str(), "%c", &cin);
    if( cin!='Y' && cin!='y' ) exit(0);

}



/*-------------------------------------- Channellist ----------------------------------------*/
/* load in channel list from input channel info */
void Channellist::Load( const std::string& chinfo )
{
    std::stringstream ss( chinfo );
    std::string channel;
    while( (ss >> channel) && (channel.at(0) != '#') )
    {
        icurrent = find(list.begin(), list.end(), channel );
        if( icurrent != list.end() )
        {
            std::cerr<<"Warning(Channellist::Load): "<<channel<<" already in the list. Will be ignored!"<<std::endl;
            continue;
        }
        list.push_back(channel);
    }
    icurrent = list.begin();

    std::cout<<"Channellist::Load: "<<list.size()<<" channels loaded ( ";
    for(unsigned int i=0; i<list.size(); i++) std::cout<<list.at(i)<<" ";
    std::cout<<")"<<std::endl;
}


/*-------------------------------------- Seedlist ----------------------------------------*/
/* load in seed list from input file and sort it by date */
static bool SameDate(const SeedInfo& a, const SeedInfo& b)
{
    return ( a.year==b.year && a.month==b.month && a.day==b.day  );
}
static bool CompareDate(const SeedInfo& a, const SeedInfo& b)
{
    return ( ( a.year<b.year ) || ( a.year==b.year && (a.month<b.month||(a.month==b.month&&a.day<b.day)) ) );
}
void Seedlist::Load( const std::string& fname )
{
    //list = new std::vector<SeedInfo>;
    std::ifstream fseed(fname);
    if( !fseed )
    {
        std::cerr<<"ERROR(Seedlist::Load): Cannot open file "<<fname<<std::endl;
        exit(-1);
    }
    std::string buff;
    SeedInfo SRtmp;
    for(; std::getline(fseed, buff);)
    {
        char stmp[buff.length()];
        if( (sscanf(buff.c_str(),"%s %d %d %d", stmp, &(SRtmp.year), &(SRtmp.month), &(SRtmp.day))) != 4 )
        {
            std::cerr<<"Warning(Seedlist::Load): format error in file "<<fname<<std::endl;
            continue;
        }
        SRtmp.name = stmp;
        list.push_back(SRtmp);
        //std::cerr<<list.back()<<std::endl;
    }
    fseed.close();
    std::cout<<"Seedlist::Load: "<<list.size()<<" seeds loaded"<<std::endl;
    /* sort the list by date */
    std::stable_sort(list.begin(), list.end(), CompareDate);
    //for(int i=0; i<list.size(); i++) std::cerr<<list.at(i)<<std::endl;
    icurrent = list.begin();
}

/* Move icurrent to (or to after) the next match of the input SeedInfo
   assuming a sorted list and relocating using binary search */
bool Seedlist::ReLocate( int year, int month, int day )
{
    SeedInfo srkey("", year, month, day);
    icurrent = std::lower_bound(list.begin(), list.end(), srkey, CompareDate );
    if( icurrent>=list.end() || icurrent<list.begin() ) return false;
    if( !SameDate(*icurrent, srkey) ) return false;
    return true;
}


/*-------------------------------------- Stationlist ----------------------------------------*/
/* load in station list from input file */
struct StaFinder
{
    StaInfo a;
    StaFinder(StaInfo b) : a(b) {}
    bool operator()(StaInfo b)
    {
        return a.name.compare(b.name)==0 && ( a.net.empty() || b.net.empty() || a.net.compare(b.net)==0 );
    }
};
bool isInt( const std::string& s ) {
   if( s.empty() ) return false;
   if( ! (isdigit(s[0]) || s[0]=='-') ) return false;
   for(int i=1; i<s.size(); i++)
      if(! isdigit(s[i])) return false;
   return true;
}
void Stationlist::Load( const std::string& fname )
{
    //list = new std::vector<StaInfo>;
    std::ifstream fsta(fname);
    if( !fsta )
    {
        std::cerr<<"ERROR(Stationlist::Load): Cannot open file "<<fname<<std::endl;
        exit(-1);
    }
    std::string buff;
    StaInfo SRtmp;
	 for(; std::getline(fsta, buff);) {
		 std::stringstream ss(buff);
		 SRtmp.CCflag = 1;
		 ss >> SRtmp.name >> SRtmp.lon >> SRtmp.lat >> SRtmp.net >> SRtmp.CCflag;
		 if( isInt(SRtmp.net) ) { SRtmp.CCflag = atoi(SRtmp.net.c_str()); SRtmp.net.clear(); }
		 icurrent = find_if(list.begin(), list.end(), StaFinder(SRtmp) );
		 if( icurrent != list.end() )
		 {
			 if( *icurrent == SRtmp )
			 {
				 std::cerr<<"Warning(Stationlist::Load):"<< SRtmp.net<<" "<<SRtmp<<" already in the list. Will be ignored!"<<std::endl;
				 continue;
			 }
			 else
			 {
				 std::cerr<<"Error(Stationlist::Load): station name confliction detected: "<<(*icurrent).net <<" "<<*icurrent<<" - "<<SRtmp.net<<" "<<SRtmp<<std::endl;
				 exit(0);
			 }
		 }
		 //std::cerr<<SRtmp.name<<" "<<SRtmp.lon<<" "<<SRtmp.lat<<" "<<SRtmp.net<<" "<<SRtmp.CCflag<<std::endl;
		 list.push_back(SRtmp);
		 //std::cerr<<list.back().name<<" "<<list.back().lon<<" "<<list.back().lat<<std::endl;
	 }
	 fsta.close();
	 std::cout<<"Stationlist::Load: "<<list.size()<<" stations loaded"<<std::endl;
	 icurrent = list.begin();
}

/* Move icurrent to the next match of the input StaInfo
	icurrent=.end() if no such match is found */
bool Stationlist::ReLocate( const std::string& staname )
{
	StaInfo srkey(staname.c_str(), 0., 0.);
	icurrent = find_if(list.begin(), list.end(), StaFinder(srkey) );
	if( icurrent>=list.end() || icurrent<list.begin() ) return false;
	return true;
}

