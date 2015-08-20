#include "SacRec.h"
#include "CCDatabase.h"
#include "CCRec.h"
#include "MyOMP.h"
#include <iostream>
#include <deque>
#include <unistd.h>

bool FileExists(const char* filename);
bool FileExists(const std::string& filename);
bool doCor(SacRec & out_CC, const SacRec & sac_am1, const SacRec & sac_ph1, const SacRec & sac_am2, const SacRec & sac_ph2,
           const int & lagtime, const int & fs, const int & mintlen, const float & tlen, const int & ftlen, std::string & fname1, std::string & fname2)   ;

std::vector < std::string > GchannelList(CCDatabase & CDB)
{
    std::istringstream ins(CDB.GetParams().chlst_info);
    // load words to this container
    std::string out;
    std::vector < std::string > channel;
    // read the words until some data in the input stream
    while (ins.good())
    {
        getline(ins, out, ' '); // tell getline to stop on ' ' character
        channel.push_back(out);
    }
    channel.pop_back();
    return channel;
}

void GStaMonList( std::vector < Station > & StaList, std::vector < std::string > & MDIR, CCDatabase & CDB)
{
    CDB.Rewind();
    int dyear = 0, dmonth, dday;
    while (1)   // Loop for station and month list
    {
        bool got;
        std::vector<DailyInfo> dinfoV;
        got = CDB.GetRec_AllCH(dinfoV);
        CDB.NextEvent();
        if (!got) break;
        bool newstationflag = true;
        for (int s = 0; s < StaList.size(); s++)
        {
            if (StaList[s].sta == dinfoV[0].staname)
            {
                newstationflag = false;
                break;
            }
        }
        if (newstationflag == true)
            StaList.push_back(Station(dinfoV[0].staname, dinfoV[0].sta.CCflag));
        bool newmonthdir=true;
        for (int m = 0; m < MDIR.size(); m++)
        {
            if (MDIR[m] == dinfoV[0].monthdir)
            {
                newmonthdir = false;
                break;
            }
        }
        if (newmonthdir == true)
            MDIR.push_back(dinfoV[0].monthdir);
    }
}

void GCCtodoList(std::vector < Station > & StaList, std::vector < std::string > & MDIR, std::vector<CC_todo> & CC_List, std::vector < std::string > & CHAN, int skipccflag)
{
    const int CHsize=CHAN.size();
    const int CCsize=CHsize*CHsize;
    for (int s1=0; s1<StaList.size(); s1++)
        for (int s2=0; s2<StaList.size(); s2++)
        {
            if (StaList[s1].sta>StaList[s2].sta || !StaList[s1].checkdoCC(StaList[s2]) ) continue;
            for (int m=0; m<MDIR.size(); m++)
            {
                bool skip_mcc=true;
                for (int f=0; f<CCsize; f++)
                {
                    int ch1=(int)(f/CHsize);
                    int ch2=f%CHsize;
                    std::string month_fout=MDIR[m]+"/COR/"+StaList[s1].sta+"/COR_" +
                                           StaList[s1].sta + "_" +CHAN[ch1] + "_" + StaList[s2].sta + "_" + CHAN[ch2] +".SAC";
                    if (!FileExists(month_fout))
                    {
                        skip_mcc=false;
                        break;
                    }
                }
                if(skip_mcc == true && skipccflag == 1)
                {
                    std::cout<<"CC File: "<<MDIR[m]<<":"<<StaList[s1].sta<<"_"<<StaList[s2].sta<<" Exist! "<<std::endl;
                    continue;
                }
                for (int d=1; d<=31; d++)
                {
                    bool all_existflag=true;
                    std::vector <std::string> infname1;
                    std::vector <std::string> infname2;
                    for (int ch=0; ch<CHsize; ch++)
                    {
                        std::string day=MDIR[m]+"."+std::to_string(d);
                        infname1.push_back(MDIR[m]+"/"+day+"/ft_"+day+"."+StaList[s1].sta+"."+CHAN[ch]+".SAC");
                        infname2.push_back(MDIR[m]+"/"+day+"/ft_"+day+"."+StaList[s2].sta+"."+CHAN[ch]+".SAC");
                        if ( !FileExists(infname1[ch]+".am") || ! FileExists(infname1[ch]+".ph") || !FileExists(infname2[ch]+".am") || !FileExists(infname2[ch]+".ph"))
                        {
                            all_existflag=false;
                            break;
                        }
                    }
                    if (all_existflag==false) continue;
                    CC_List.push_back( CC_todo(infname1,infname2, MDIR[m], StaList[s1].sta, StaList[s2].sta, d) );
                }
            }
        }
}

void CCList2CC( std::vector< CC_todo > & CC_List, std::vector < std::string > & CHAN, CCDatabase & CDB)
{
    const int CHsize = CHAN.size();
    const int CCsize = CHsize*CHsize;
    int fskipcc=CDB.GetParams().fskipcrco;
    int ftlen=CDB.GetParams().ftlen;
    float tlen=CDB.GetParams().tlen;
    int mintlen=CDB.GetParams().mintlen;
    int SIZE = CC_List.size();
    if (fskipcc !=2 && SIZE!=0)
    {
        std::cout<<"ATTENTION: Start Cross Correlation for "<<SIZE<<" pairs!"<<std::endl;
        int i=-1;
        int si=-1;
        const int fs = CDB.GetParams().sps;
        const int LagTime = CDB.GetParams().lagtime;
        const int Coutflag=CDB.GetParams().CorOutflag;
        const int checkprec=CDB.GetParams().fprcs;
        std::vector < outCCRec > CCRec_list; // 1.06--> 1.07
        CCRec_list.reserve(50);
        #pragma omp parallel private ( i ) shared( CCRec_list, CC_List, si)
        {
            #pragma omp sections nowait
            {
                #pragma omp section
                {
                    int iconsume=0;
                    bool newflag=true;
                    std::vector < SacRec > m_sacV(CCsize);
                    std::vector < std::string> foutnameM(CCsize);
                    std::string MDIR=CC_List[0].monthdir; // Initiation
                    std::string STA1=CC_List[0].sta1;
                    std::string STA2=CC_List[0].sta2;
                    while(1)
                    {
                        if (iconsume>=SIZE) break;
                        bool skipcon=false;
                        #pragma omp critical (CCRecIO)
                        {
                            if (CCRec_list.empty()) skipcon=true;
                        }
                        if (skipcon==true)
                        {
                            sleep(0.1);
                            continue;
                        }
                        #pragma omp critical (CCRecIO)
                        {
                            if (iconsume+1!=CCRec_list[0].cclistNo) skipcon=true;
                        }
                        if (skipcon==true)
                        {
                            sleep(0.1);
                            continue;
                        }
                        iconsume++;
                        if (MDIR!=CCRec_list[0].monthdir || STA1!=CCRec_list[0].sta1 || STA2!=CCRec_list[0].sta2)
                        {
                            newflag=true; // Update!
                            MDIR=CCRec_list[0].monthdir;
                            STA1=CCRec_list[0].sta1;
                            STA2=CCRec_list[0].sta2;
                            for (int j=0; j<CCsize; j++)
                            {
                                auto & m_sac=m_sacV[j];
                                if (! m_sac.sig) continue;
                                m_sac.Write(foutnameM[j]);
                            }
                        }
                        if (newflag==true)
                        {
                            #pragma omp critical (CCRecIO)
                            {
                                MKDirs((CCRec_list[0].COR).c_str());
                                if (CCRec_list[0].CCRec.size()==CCsize)
                                    for (int j=0; j<CCsize; j++)
                                    {
                                        if (! CCRec_list[0].CCRec[j].sig ) break;
                                        m_sacV[j] = std::move( CCRec_list[0].CCRec[j]) ;
                                        foutnameM[j]=CCRec_list[0].outfname[j];
                                        newflag=false;
                                    }
                            }
                        }
                        else
                        {
                            if ( CCRec_list[0].CCRec.size() == CCsize )
                                for (int j=0; j<CCsize; j++) // Monthly stacking
                                {
                                    auto & inSac=CCRec_list[0].CCRec[j];
                                    if (! m_sacV[j].sig || ! inSac.sig) break;
                                    m_sacV[j].shd.user0+=1;
                                    for (int k=0; k<inSac.shd.npts; k++)
                                        m_sacV[j].sig[k]+=inSac.sig[k];
                                }
                        }
                        #pragma omp critical (CCRecIO)
                        {
                             CCRec_list.erase(CCRec_list.begin());
                        }
                    }
                    for (int j=0; j<CCsize; j++)
                    {
                        auto & m_sac=m_sacV[j];
                        if (! m_sac.sig) continue;
                        m_sac.Write(foutnameM[j]);
                    }
                }
            }
            while (1)
            {
                #pragma omp critical (ListNoinc)
                {
                    si++;
                    i=si;
                }
                if(i>=SIZE) break;
                std::vector < SacRec > d_sacV(CCsize);
                std::vector < std::string> foutnameD(CCsize);
                std::vector <std::string> m_foutname(CCsize);
                std::vector < SacRec > sac_amV1(CHsize), sac_phV1(CHsize), sac_amV2(CHsize), sac_phV2(CHsize);
                std::string cor_dir=CC_List[i].monthdir + "/COR/" + CC_List[i].sta1;
                bool skipccflag=false;
                for (int ch=0; ch<CHsize; ch++)
                {
                    sac_amV1[ch].Load(CC_List[i].infname1[ch]+".am");
                    sac_phV1[ch].Load(CC_List[i].infname1[ch]+".ph");
                    sac_amV2[ch].Load(CC_List[i].infname2[ch]+".am");
                    sac_phV2[ch].Load(CC_List[i].infname2[ch]+".ph");
                }
                for (int j=0; j<CCsize; j++)
                {
                    int ch1=(int)(j/CHsize); // 1.05-->1.06
                    int ch2=j%CHsize;
                    foutnameD[j]= CC_List[i].monthdir+ "/COR_D/" + CC_List[i].sta1 + "/COR_" +
                                  CC_List[i].sta1 + "_" +CHAN[ch1] + "_" + CC_List[i].sta2 + "_" + CHAN[ch2]  +
                                  "_" + std::to_string(CC_List[i].day) + ".SAC";
                    m_foutname[j]=CC_List[i].monthdir+ "/COR/" + CC_List[i].sta1 + "/COR_" +
                                  CC_List[i].sta1 + "_" +CHAN[ch1] + "_" + CC_List[i].sta2 + "_" + CHAN[ch2] +".SAC";
                    auto & d_sac=d_sacV[j];
                    auto & foutname=foutnameD[j];
                    // Start to do CC
                    d_sac.LoadHD(CC_List[i].infname1[ch1]+".am");
                    std::string frec1=CC_List[i].infname1[ch1]+"_rec";
                    std::string frec2=CC_List[i].infname1[ch2]+"_rec";
                    if (doCor(d_sac, sac_amV1[ch1], sac_phV1[ch1], sac_amV2[ch2], sac_phV2[ch2], LagTime, fs, mintlen, tlen, ftlen, frec1, frec2))
                        std::cout << " THREAD: ("<<omp_get_thread_num()<<") ---- Do CC: " << CC_List[i].sta1 << "_" << CHAN[ch1]  +
                                  " with " + CC_List[i].sta2 + "_" + CHAN[ch2] << " Date: " <<
                                  CC_List[i].monthdir << "." <<CC_List[i].day <<" ---- List No. "<<i+1<< std::endl;
                    else
                    {
                        std::cout << " THREAD: ("<<omp_get_thread_num()<<") ---- Error am&ph File CC: " << CC_List[i].sta1 << "_" << CHAN[ch1]  +
                                  " with " + CC_List[i].sta2 + "_" + CHAN[ch2] << " Date: " <<
                                  CC_List[i].monthdir << "." <<CC_List[i].day <<" ---- List No. "<<i+1<< std::endl;
                        skipccflag=true;
                        break;
                    }
                    if (checkprec==1) // NEED TEST!
                    {
                        if (d_sac.CheckPrecNoise()) skipccflag=true;
                        break;
                    }
                    if (Coutflag != 0 )
                    {
                        MKDirs((CC_List[i].monthdir + "/COR_D/" + CC_List[i].sta1).c_str());
                        d_sac.Write(foutname); // Save daily CC data
                    }
                }
                if (skipccflag==true)
                {
                    for (int j=0; j<CCsize; j++)
                        if (d_sacV[j].sig) d_sacV[j].sig.reset();
                }
                if(Coutflag!=1)
                {
                    #pragma omp critical (CCRecIO)
                    {
                        int N=i+1;
                        if (CCRec_list.size()==0)
                            CCRec_list.push_back(outCCRec(d_sacV, CCsize, CC_List[i].day, N, CC_List[i].monthdir,CC_List[i].sta1, CC_List[i].sta2, m_foutname, cor_dir));
                        else
                        {
                            std::vector < outCCRec >:: iterator InsertStart=CCRec_list.end()-1; // 1.06--> 1.07
                            while ((*InsertStart).cclistNo>N && InsertStart!=CCRec_list.begin())
                            {
                                InsertStart--;
                            }
                            if(InsertStart==CCRec_list.begin() && (*InsertStart).cclistNo>N)
                                CCRec_list.insert(CCRec_list.begin(),outCCRec(d_sacV,  // 1.06--> 1.07
                                                                              CCsize, CC_List[i].day, N, CC_List[i].monthdir, CC_List[i].sta1, CC_List[i].sta2, m_foutname, cor_dir));
                            else
                                CCRec_list.insert(++InsertStart, outCCRec(d_sacV, CCsize, CC_List[i].day, N, CC_List[i].monthdir, CC_List[i].sta1, CC_List[i].sta2, m_foutname, cor_dir));
                        }
                    }
                }
            }
        }
    }
}

