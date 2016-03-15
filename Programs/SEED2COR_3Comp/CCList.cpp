#include "SacRec.h"
#include "CCDatabase.h"
#include "CCRec.h"
#include "MyOMP.h"
#include <iostream>
#include <deque>
#include <unistd.h>
#include "StaPair.h"

bool FileExists(const char* filename);
bool FileExists(const std::string& filename);
void GCCtodoList( std::vector<StaPair> & CC_List, std::vector < StaInfo > & StaList, std::vector < std::string > & MDIR, std::vector < std::vector < std::string > > & CHAll, int skipccflag)
{
    const int CHsize=CHAll.size();
    const int CCsize=CHsize*CHsize;
    for (int s1=0; s1<StaList.size(); s1++)
        for (int s2=0; s2<StaList.size(); s2++)
        {
            if (StaList[s1].name>StaList[s2].name || !StaList[s1].checkdoCC(StaList[s2]) ) continue;
            StaPair tempSPair=StaPair( StaList[s1], StaList[s2] );
            tempSPair.GetDistAzBaz();
            for (int m=0; m<MDIR.size(); m++)
            {
                bool skip_mcc=true;
                for (int f=0; f<CCsize; f++)
                {
                    int ENo=0;
                    int ch1=(int)(f/CHsize);
                    int ch2=f%CHsize;
                    for (int index1=0; index1<CHAll[ch1].size(); index1++)
                        for (int index2=0; index2<CHAll[ch2].size(); index2++)
                        {
                            std::string month_fout=MDIR[m]+"/COR/"+StaList[s1].name+"/COR_" +
                                                   StaList[s1].name + "_" +CHAll[ch1][index1] + "_" + StaList[s2].name + "_" + CHAll[ch2][index2] +".SAC";
                            if (!FileExists(month_fout))
                            {
                                skip_mcc=false;
                                break;
                            }
                            ENo++;
                        }
                    if (ENo>1)
                        std::cout<<ENo<<" of CCs exist for:"<< StaList[s1].name<<"_"<<StaList[s2].name<<std::endl;
                }

                if(skip_mcc == true && skipccflag == 1)
                {
                    std::cout<<"CC File: "<<MDIR[m]<<":"<<StaList[s1].name<<"_"<<StaList[s2].name<<" Exist! "<<std::endl;
                    continue;
                }

                for (int d=1; d<=31; d++)
                {
                    std::vector < std::string > tempCH1, tempCH2;
                    std::vector <std::string> infname1;
                    std::vector <std::string> infname2;
                    std::string day=MDIR[m]+"."+std::to_string(d);
                    for (int ch=0; ch<CHsize; ch++)
                    {
                        for (int i=0; i<CHAll[ch].size(); i++)
                        {
                            std::string infname1_pre=MDIR[m]+"/"+day+"/ft_"+day+"."+StaList[s1].name+"."+CHAll[ch][i]+".SAC";
                            if ( FileExists(infname1_pre+".am") && FileExists(infname1_pre+".ph") )
                            {
                                tempCH1.push_back(CHAll[ch][i]);
                                break;
                            }
                        }
                        for (int i=0; i<CHAll[ch].size(); i++)
                        {
                            std::string infname2_pre=MDIR[m]+"/"+day+"/ft_"+day+"."+StaList[s2].name+"."+CHAll[ch][i]+".SAC";
                            if ( FileExists(infname2_pre+".am") && FileExists(infname2_pre+".ph") )
                            {
                                tempCH2.push_back(CHAll[ch][i]);
                                break;
                            }
                        }
                    }
                    if (tempCH1.size()!=0 && tempCH2.size()!=0)   ///
                    {
                        tempSPair.SetDate(MDIR[m], d);
                        tempSPair.SetChanLst(tempCH1, tempCH2);
                        CC_List.push_back( tempSPair );
                    }
                }
            }
        }
}

void CCList2CC( std::vector<StaPair> & CC_List, std::vector < std::vector < std::string > > ChanAll, CCDatabase & CDB)
{
    const int CHsize = ChanAll.size();
    const int CCsize = CHsize*CHsize;
    int fskipcc=CDB.GetParams().fskipcrco;
    int SIZE = CC_List.size();
    if (fskipcc !=2 && SIZE!=0)
    {
        std::cout<<"ATTENTION: Start Cross Correlation for "<<SIZE<<" pairs!"<<std::endl;
        int i=-1;
        int si=-1;
        const int fs = CDB.GetParams().sps;
        const int Coutflag=CDB.GetParams().CorOutflag;
        std::vector < outCCRec > CCRec_list; // 1.06--> 1.07
        std::vector < std::string > doDateLst, errorLst;
        std::string MDIR, STA1, STA2; // Initiation ###
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
                    /// std::string MDIR, STA1, STA2;
                    std::vector < std::string >  chanL1(CHsize), chanL2(CHsize);
                    bool initial_flag=false;
                    while(1)
                    {
                        if (iconsume>=SIZE || Coutflag==1) break;
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
                        if (initial_flag ==false && CCRec_list[0].sacRecSize==CCsize )
                        {
                            MDIR=CC_List[0].month; // Initiation
                            STA1=CC_List[0].sta1.name;
                            STA2=CC_List[0].sta2.name;
                            for (int i=0; i<CC_List[0].chanL1.size(); i++)
                                chanL1[i]=CC_List[0].chanL1[i];
                            for (int i=0; i<CC_List[0].chanL2.size(); i++)
                                chanL2[i]=CC_List[0].chanL2[i];
                            initial_flag=true;
                        }
                        iconsume++;
                        if (initial_flag ==true)
                        {
                            for (int i=0; i<CC_List[0].chanL1.size(); i++)
                            {
                                if (chanL1[i]!=CC_List[0].chanL1[i])
                                {
                                    newflag=true;
                                    break;
                                }
                            }

                            for (int i=0; i<CC_List[0].chanL2.size(); i++)
                            {
                                if (chanL2[i]!=CC_List[0].chanL2[i])
                                {
                                    newflag=true;
                                    break;
                                }
                            }

                            if (MDIR!=CCRec_list[0].monthdir || STA1!=CCRec_list[0].sta1 || STA2!=CCRec_list[0].sta2)
                            {
                                newflag=true; // Update!
                                if (STA1!=CCRec_list[0].sta1 || STA2!=CCRec_list[0].sta2)
                                {
                                    #pragma omp critical (PrintCC)
                                    {
                                        bool outflag=false;
                                        if (doDateLst.size()!=0)
                                        {
                                            outflag=true;
                                            std::cout<<"Do CC for: "<<STA1<<"_"<<STA2<<":"<<std::endl;
                                            for (int i=0; i<doDateLst.size(); i++)
                                                std::cout<<" "<<doDateLst[i]<<" ";
                                            std::cout<<" "<<std::endl;
                                            doDateLst.clear();
                                        }
                                        if (errorLst.size()!=0)
                                        {
                                            outflag=true;
                                            std::cout<<"Error CC for: "<<STA1<<"_"<<STA2<<":"<<std::endl;
                                            for (int i=0; i<errorLst.size(); i++)
                                                std::cout<<" "<<errorLst[i]<<" ";
                                            std::cout<<" "<<std::endl;
                                            errorLst.clear();
                                        }
                                        if (outflag==true)
                                            std::cout<<"================================================================="<<std::endl;
                                    }
                                }
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
                std::vector < SacRec > d_sacV;
                std::vector < std::string> foutnameD;
                std::vector <std::string> m_foutname;
                CC_List[i].staPair2CC( doDateLst, errorLst, d_sacV, foutnameD, m_foutname, CDB, i);
                int d_CCsize=d_sacV.size();
                std::string cor_dir=CC_List[i].month + "/COR/" + CC_List[i].sta1.name;
                if( Coutflag!=1 )
                {
                    #pragma omp critical (CCRecIO)
                    {
                        int N=i+1;
                        if (CCRec_list.size()==0)
                            CCRec_list.push_back(outCCRec(d_sacV, d_CCsize, CC_List[i].day, N, CC_List[i].month, CC_List[i].sta1.name, CC_List[i].chanL1,
                            CC_List[i].sta2.name, CC_List[i].chanL2, m_foutname, cor_dir));
                        else
                        {
                            std::vector < outCCRec >:: iterator InsertStart=CCRec_list.end()-1; // 1.06--> 1.07
                            while ((*InsertStart).cclistNo>N && InsertStart!=CCRec_list.begin())
                            {
                                InsertStart--;
                            }
                            if(InsertStart==CCRec_list.begin() && (*InsertStart).cclistNo>N)
                                CCRec_list.insert(CCRec_list.begin(),outCCRec(d_sacV, d_CCsize, CC_List[i].day, N, CC_List[i].month, CC_List[i].sta1.name, CC_List[i].chanL1,
                                CC_List[i].sta2.name, CC_List[i].chanL2, m_foutname, cor_dir));
                            else
                                CCRec_list.insert(++InsertStart, outCCRec(d_sacV, d_CCsize, CC_List[i].day, N, CC_List[i].month, CC_List[i].sta1.name, CC_List[i].chanL1,
                                CC_List[i].sta2.name, CC_List[i].chanL2, m_foutname, cor_dir));
                        }
                    }
                }
            }

            #pragma omp critical (PrintCC)
            {
                bool outflag=false;
                if (doDateLst.size()!=0)
                {
                    outflag=true;
                    std::cout<<"Do CC for: "<<STA1<<"_"<<STA2<<":"<<std::endl;
                    for (int i=0; i<doDateLst.size(); i++)
                        std::cout<<" "<<doDateLst[i]<<" ";
                    std::cout<<" "<<std::endl;
                    doDateLst.clear();
                }
                if (errorLst.size()!=0)
                {
                    outflag=true;
                    std::cout<<"Error CC for: "<<STA1<<"_"<<STA2<<":";
                    for (int i=0; i<errorLst.size(); i++)
                        std::cout<<" "<<errorLst[i]<<" ";
                    std::cout<<" "<<std::endl;
                    errorLst.clear();
                }
                if (outflag==true)
                    std::cout<<"================================================================="<<std::endl;
            }
        }
    }
}

