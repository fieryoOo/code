#include <iostream>
#include <vector>
#include <string>
#include "StaPair.h"
#include "MyLogger.h"
#include "SacRec.h"
#include "CCDatabase.h"
#include "CCRec.h"

extern MyLogger logger;

double vincenty_earth_dist( const double lat1, const double lon1,
                            const double lat2, const double lon2,
                            const double flat, double *azimuth,
                            double *back_azimuth);
bool doCor(SacRec & out_CC, const SacRec & sac_am1, const SacRec & sac_ph1, const SacRec & sac_am2, const SacRec & sac_ph2,
           const int & lagtime, const int & fs, const int & mintlen, const float & tlen, const int & ftlen, std::string & fname1, std::string & fname2,
           const double dist, const double azimuth, const double back_azimuth);  // Note that lagtime and fs are integer!!!
StaPair::StaPair() {}
StaPair::~StaPair() {}
StaPair::StaPair(const StaInfo Sta1, const StaInfo Sta2)
{
    sta1=Sta1;
    sta2=Sta2;
    dist=-1;
}

bool StaPair::GetDistAzBaz()
{
    const double lat1=sta1.lat;
    const double lon1=sta1.lon;
    const double lat2=sta2.lat;
    const double lon2=sta2.lon;
    if (abs(lat1-lat2) <0.001 && abs(lon1-lon2) <0.001) return false;
    const double flattening = 1. / 298.257223563;
    const double PI = 3.14159265358979323846;
    double azimuth, DIST, back_azimuth,evla,evlo, stla, stlo;
    evla=lat1*PI/180.0;
    evlo=lon1*PI/180.0;
    stla=lat2*PI/180.0;
    stlo=lon2*PI/180.0;
    DIST = vincenty_earth_dist( evla, evlo, stla, stlo, flattening, &azimuth,
                                &back_azimuth);
    dist=DIST;
    az=azimuth;
    baz=back_azimuth;
    return true;
}

void StaPair::SetDate( const std::string &Month, const int & Day )
{
    month = Month;
    day=Day;
    return;
}

void StaPair::SetChanLst(const std::vector < std::string >  ChanL1, const std::vector < std::string > ChanL2)
{
    chanL1.clear();
    chanL2.clear();
    for (int i=0; i<ChanL1.size(); i++ )
    {
        chanL1.push_back(ChanL1[i]);
    }
    for (int i=0; i<ChanL2.size(); i++ )
    {
        chanL2.push_back(ChanL2[i]);
    }
}


void StaPair::staPair2CC(std::vector < SacRec > & D_sacV, std::vector < std::string> & FoutnameD,
                         std::vector <std::string> & M_foutname, CCDatabase & cdb, int LNo )
{
    const int CHsize1= chanL1.size();
    const int CHsize2= chanL2.size();
    int fskipcc=cdb.GetParams().fskipcrco;
    int ftlen=cdb.GetParams().ftlen;
    float tlen=cdb.GetParams().tlen;
    int mintlen=cdb.GetParams().mintlen;
    const int fs = cdb.GetParams().sps;
    const int LagTime = cdb.GetParams().lagtime;
    const int Coutflag=cdb.GetParams().CorOutflag;
    const int checkprec=cdb.GetParams().fprcs;
    std::vector < SacRec > sac_amV1(CHsize1), sac_phV1(CHsize1), sac_amV2(CHsize2), sac_phV2(CHsize2); ///

    std::string MDAY=month+"."+std::to_string(day);
    for (int ch=0; ch<CHsize1; ch++)
    {
        std::string infname1=month+"/"+MDAY+"/ft_"+MDAY+"."+sta1.name+"."+chanL1[ch]+".SAC";
        sac_amV1[ch].Load(infname1+".am");
        sac_phV1[ch].Load(infname1+".ph");
    }
    for (int ch=0; ch<CHsize2; ch++)
    {
        std::string infname2=month+"/"+MDAY+"/ft_"+MDAY+"."+sta2.name+"."+chanL2[ch]+".SAC";
        sac_amV2[ch].Load(infname2+".am");
        sac_phV2[ch].Load(infname2+".ph");
    }
    bool skipccflag=false;
    for (int ch1=0; ch1<CHsize1; ch1++)
        for (int ch2=0; ch2<CHsize2; ch2++)
        {
            std::string foutname=month+ "/COR_D/" + sta1.name + "/COR_" +sta1.name + "_" +chanL1[ch1] + "_"
                + sta2.name + "_" + chanL2[ch2]  +"_" + std::to_string(day) + ".SAC";
            FoutnameD.push_back(foutname);
            M_foutname.push_back(month+ "/COR/" + sta1.name + "/COR_" +
                                 sta1.name + "_" +chanL1[ch1] + "_" + sta2.name + "_" + chanL2[ch2]  + ".SAC");
            SacRec d_sac;
            // Start to do CC
            std::string infname1=month+"/"+MDAY+"/ft_"+MDAY+"."+sta1.name+"."+chanL1[ch1]+".SAC";
            std::string infname2=month+"/"+MDAY+"/ft_"+MDAY+"."+sta2.name+"."+chanL2[ch2]+".SAC";
            d_sac.LoadHD(infname1+".am");
            std::string frec1=infname1+"_rec";
            std::string frec2=infname2+"_rec";
            if (doCor(d_sac, sac_amV1[ch1], sac_phV1[ch1], sac_amV2[ch2], sac_phV2[ch2], LagTime, fs, mintlen, tlen, ftlen, frec1, frec2, dist, az, baz))
                std::cout << " THREAD: ("<<omp_get_thread_num()<<") ---- Do CC: " << sta1.name << "_" << chanL1[ch1]  +
                          " with " + sta2.name + "_" + chanL2[ch2] << " Date: " <<
                          month<< "." <<day <<" ---- List No. "<<LNo+1<< std::endl;
            else
            {
                std::cout << " THREAD: ("<<omp_get_thread_num()<<") ---- Error am&ph File CC: " << sta1.name << "_" << chanL1[ch1]  +
                          " with " + sta2.name + "_" + chanL2[ch2] << " Date: " <<
                          month<< "." <<day <<" ---- List No. "<<LNo+1<< std::endl;
                skipccflag=true;
                break;
            }
            if (checkprec==1) // NEED TEST!
            {
                if (d_sac.CheckPrecNoise()) skipccflag=true;
                break;
            }  ///
            if (Coutflag != 0 )
            {
                MKDirs((month+ "/COR_D/" + sta1.name).c_str());
                d_sac.Write(foutname); // Save daily CC data
            }
            if (skipccflag==false) D_sacV.push_back(d_sac);
        }
}
