#include "SacRec.h"
#include "SysTools.h"
#include "CCRec.h"
#include "MyLogger.h"
#include "MyOMP.h"
#include "CCDatabase.h"
#include <fftw3.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <limits>
#include <chrono>
#include <random>
#include <algorithm>
#include <sys/stat.h>
#include <complex.h>
#include <deque>
extern MyLogger logger;

double vincenty_earth_dist( const double lat1, const double lon1,
                            const double lat2, const double lon2,
                            const double flat, double *azimuth,
                            double *back_azimuth);
bool CalcRecCor( std::string & fname1, std::string fname2, float *cor_rec, int lagn, int mintlen, float tlen, float fs);
//#include <pthread.h>

bool FileExists(const char* filename)
{
    struct stat info;
    int ret = -1;

    //get the file attributes
    ret = stat(filename, &info);
    if(ret == 0)
    {
        //stat() is able to get the file attributes,
        //so the file obviously exists
        return true;
    }
    else
    {
        //stat() is not able to get the file attributes,
        //so the file obviously does not exist or
        //more capabilities is required
        return false;
    }
}
bool FileExists(const std::string& filename)
{
    return FileExists(filename.c_str());
}
void Rotation_daily(std::deque < SacRec > & SACV, std::vector < std::string > Rout)
{
    // Component Correction for daily data. Need TEST!!!
    if (SACV.size()!=Rout.size() || abs ( (SACV[0].shd.cmpaz-90)- SACV[1].shd.cmpaz)>0.1 )
    {
        std::cout<<"WARNING: Incompatible daily rotation input!!"<<std::endl;
        return;
    }

    std::cout<< "ATTENTION: Daily component correction: "<<Rout[0]<< " and "<<Rout[1]<< std::endl;
    const float PI=3.14159265358979323846;
    float del=SACV[1].shd.cmpaz/180*PI;
    auto & Corr_E=SACV[0];
    auto & Corr_N=SACV[1]; // Corrected E and N component
    float tempN,tempE;
    for (int i=0; i<SACV[0].shd.npts; i++)
    {
        tempN=Corr_N.sig[i];
        tempE=Corr_E.sig[i];
        Corr_N.sig[i]=tempN*cos(del)-tempE*sin(del);
        Corr_E.sig[i]=tempN*sin(del)+tempE*cos(del);
    }
    Corr_E.shd.cmpaz=90.0;
    Corr_N.shd.cmpaz=0.0;
    Corr_E.Write(Rout[0]);
    Corr_N.Write(Rout[1]);
}
void stackCCdata(CC_output & CCstack, int checkprec)
{
    std::string stackdir;
    SacRec StackedData;
    StackedData.Load(CCstack.outfname[0]);
    //StackedData.shd.nzjday=1;
    //logger.Hold( INFO, "Stack all for Sta & Chan pairs: " + CCstack.sta1 + "_" + CCstack.chan1 + "_" + CCstack.sta2 + "_" + CCstack.chan2 , FuncName );
    std::cout << "Stack all for Sta&chan pairs: " << CCstack.sta1 << "_" << CCstack.chan1 << "_" << CCstack.sta2 << "_" << CCstack.chan2 << std::endl;
    stackdir = "COR/" + CCstack.sta1;
    StackedData.shd.nzjday = 1;
    MKDirs((stackdir).c_str());
    for (int j = 1; j < CCstack.outfname.size(); j++)
    {
        SacRec singleData;
        singleData.Load(CCstack.outfname[j]);
        if (checkprec==1)
            if(singleData.CheckPrecNoise()) continue;
        if (singleData.shd.npts != StackedData.shd.npts)
        {
            logger.Hold( WARNING, "Incompatible npts for file: " + CCstack.outfname[j]  , FuncName );
            //std::cout << "Warning: incompatible npts for file: " << CCstack.outfname[j] << std::endl;
            continue;
        }
        logger.Hold( INFO, " FILE: " + CCstack.outfname[j]  , FuncName );
        // std::cout << " FILE: " << CCstack.outfname[j] << std::endl;
        StackedData.shd.user0 += singleData.shd.user0;
        for (int k = 0; k < StackedData.shd.npts; k++)
            StackedData.sig[k] += singleData.sig[k];
    }
    StackedData.Write(stackdir + "/COR_" + CCstack.sta1 + "_" + CCstack.chan1 + "_" + CCstack.sta2 + "_" + CCstack.chan2 + ".SAC");
}

CC_todo::CC_todo(const std::vector <std::string> & IN1, const std::vector <std::string> &IN2, const std::string & MDIR, const std::string & STA1,
                 const std::string & STA2,  const int & DD)
{
    infname1=IN1;
    infname2=IN2;
    monthdir=MDIR;
    sta1=STA1;
    sta2=STA2;
    day=DD;

}
CC_todo::CC_todo( ) {}
CC_todo::~CC_todo( ) {}

bool doCor(SacRec & out_CC, const SacRec & sac_am1, const SacRec & sac_ph1, const SacRec & sac_am2, const SacRec & sac_ph2,
           const int & lagtime, const int & fs, const int & mintlen, const float & tlen, const int & ftlen, std::string & fname1, std::string & fname2)  // Note that lagtime and fs are integer!!!
{
    // xCorr function:
    // Definition: (sac1 x sac2)[n] = SIGMA( sac1[m] *sac2[m+n] )
    int N=sac_am1.shd.npts;
    int Ns=2*N-1; ///??? or 2*N-1???
    std::complex<double> temp1,temp2;
    fftw_plan px;
    fftw_complex * CCsp = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Ns);
    fftw_complex * out_shifted = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Ns);
    #pragma omp critical
    {
        px = fftw_plan_dft_1d(Ns, CCsp, out_shifted, FFTW_BACKWARD, FFTW_ESTIMATE);
    }
    if(isnan(sac_ph1.sig[0])||isnan(sac_ph2.sig[0])||isnan(sac_am1.sig[0])||isnan(sac_am2.sig[0]))
    {
        //logger.Hold( WARNING, "Skip CC due to NaN!" , FuncName );
        std::cout<<"Warning: Skip CC due to NaN!"<<std::endl;
        return false;
    }
    // set value for the first point of CCsp
    temp1=std::complex<double>(sac_am1.sig[0]*cos(sac_ph1.sig[0]),sac_am1.sig[0]*sin(sac_ph1.sig[0]));
    temp2=std::complex<double>(sac_am2.sig[0]*cos(sac_ph2.sig[0]),sac_am2.sig[0]*sin(sac_ph2.sig[0]));
    CCsp[0][0]=(temp2*conj(temp1)).real();
    CCsp[0][1]=(temp2*conj(temp1)).imag();
    // In frequency domain: conj(sac1) * sac2
    for(int i=1; i<N; i++ )
    {
        temp1=std::complex<double>(sac_am1.sig[i]*cos(sac_ph1.sig[i]),sac_am1.sig[i]*sin(sac_ph1.sig[i]));
        temp2=std::complex<double>(sac_am2.sig[i]*cos(sac_ph2.sig[i]),sac_am2.sig[i]*sin(sac_ph2.sig[i]));
        CCsp[i][0]=(temp2*conj(temp1)).real();
        CCsp[i][1]=(temp2*conj(temp1)).imag();
        CCsp[i+N-1][0]=0; // zeropadding
        CCsp[i+N-1][1]=0;
        if ( sac_am1.sig[i] > 10e20 || abs(sac_ph1.sig[i])> 10e20 || sac_am2.sig[i]>10e20 || abs(sac_ph2.sig[i])>10e20)
            return false;
    }
    fftw_execute(px);
    #pragma omp critical
    {
        fftw_destroy_plan(px);
    }
    fftw_free(CCsp);
    float *seis_out;
    seis_out= (float *)malloc(sizeof(float) *Ns);  /// ver. 1.00 --> ver. 1.01
    for (int j=0; j<Ns; j++)
        seis_out[j]=2*out_shifted[j][0]; // scaling, Need modification to consider length of record, holes,etc.
    fftw_free(out_shifted);
    int lagn=(int)floor(lagtime*fs +0.5);
    if (lagn>Ns)
    {
        lagn=Ns;
        logger.Hold( WARNING, "Lagtime overflow!" , FuncName );
        //std::cout<<"Warning: lagtime overflow!" <<std::endl;
    }
    float depmin=10000;
    float depmax=-10000;

    float *sigtmp = new float[2*lagn+1];
    float sums;
    out_CC.sig=std::unique_ptr<float[]>(sigtmp);
    if (ftlen==0)
    {
        out_CC.sig[lagn]=seis_out[0]/Ns;
        sums=seis_out[0]/Ns;
        depmax=std::max(depmax,out_CC.sig[lagn]);
        depmin=std::min(depmin,out_CC.sig[lagn]);
        for(int k=1; k<(lagn+1); k++)
        {
            out_CC.sig[lagn-k]=seis_out[k]/Ns; // lagn-1 ~ 0
            out_CC.sig[lagn+k]=seis_out[Ns-k]/Ns; // lagn +1 ~ 2*lagn
            sums=sums+seis_out[k]/Ns+seis_out[Ns-k]/Ns;
            depmax=std::max(depmax,out_CC.sig[lagn-k]);
            depmax=std::max(depmax,out_CC.sig[lagn+k]);
            depmin=std::min(depmin,out_CC.sig[lagn-k]);
            depmin=std::min(depmin,out_CC.sig[lagn+k]);
        }
    }
    else
    {
        float *cor_rec;
        cor_rec= (float *)malloc(sizeof(float) *(2*lagn+1));
        if ( ! CalcRecCor( fname1,fname2, cor_rec, lagn, mintlen, tlen, fs))
            return false;
        out_CC.sig[lagn]=seis_out[0]/cor_rec[lagn];
        sums=out_CC.sig[lagn];
        depmax=std::max(depmax,out_CC.sig[lagn]);
        depmin=std::min(depmin,out_CC.sig[lagn]);
        for(int k=1; k<(lagn+1); k++)
        {
            out_CC.sig[lagn-k]=seis_out[k]/cor_rec[lagn-k]; // lagn-1 ~ 0
            out_CC.sig[lagn+k]=seis_out[Ns-k]/cor_rec[lagn+k]; // lagn +1 ~ 2*lagn
            sums=sums+out_CC.sig[lagn-k]+out_CC.sig[lagn+k];
            depmax=std::max(depmax,out_CC.sig[lagn-k]);
            depmax=std::max(depmax,out_CC.sig[lagn+k]);
            depmin=std::min(depmin,out_CC.sig[lagn-k]);
            depmin=std::min(depmin,out_CC.sig[lagn+k]);
        }
        free(cor_rec);
    }
    free(seis_out);
    /*---- set header value----*/
    out_CC.shd.npts=2*lagn+1;
    out_CC.shd.b=-(float)lagn/(float)fs;
    out_CC.shd.e= (float)lagn/(float)fs;
    out_CC.shd.iftype=1;
    out_CC.shd.leven=1;
    out_CC.shd.delta=1/(float)fs;
    out_CC.shd.idep=1;
    out_CC.shd.depmax=depmax;
    out_CC.shd.depmin=depmin;
    out_CC.shd.depmen=sums/(float)(2*lagn+1);
    out_CC.shd.nzhour=0;
    out_CC.shd.nzmin=0;
    out_CC.shd.nzsec=0;
    sscanf(sac_am1.shd.kcmpnm,"%s",out_CC.shd.kcmpnm); // kcmpnm is not just "BHE", it is "BHE     DK      -12345  -12345  <B0>"
    char temp[3];
    sscanf(sac_am2.shd.kcmpnm,"%s",temp);
    strcat(out_CC.shd.kcmpnm,temp);
    sscanf(sac_am2.shd.kstnm,"%s",out_CC.shd.kstnm);
    sscanf(sac_am1.shd.kstnm,"%s",out_CC.shd.kevnm);
    sscanf(sac_am2.shd.knetwk,"%s",out_CC.shd.knetwk);
    out_CC.shd.stlo=sac_am2.shd.stlo;
    out_CC.shd.stla=sac_am2.shd.stla;
    out_CC.shd.evlo=sac_am1.shd.stlo;
    out_CC.shd.evla=sac_am1.shd.stla;
    out_CC.shd.cmpinc=0.0; /// TODO
    out_CC.shd.stel=sac_am2.shd.stel;
    out_CC.shd.stdp=sac_am2.shd.stdp;
    out_CC.shd.lpspol=1;
    out_CC.shd.lcalda=1;
    //out_CC.shd.lovrok=1;
    // compute dist, az, baz using vincenty formula
    const double flattening = 1. / 298.257223563;
    const double semimajor = 6378.137;
    const double PI = 3.14159265358979323846;
    double azimuth, dist, back_azimuth,evla,evlo, stla, stlo;
    evla=out_CC.shd.evla*PI/180.0;
    evlo=out_CC.shd.evlo*PI/180.0;
    stla=out_CC.shd.stla*PI/180.0;
    stlo=out_CC.shd.stlo*PI/180.0;
    dist = vincenty_earth_dist( evla, evlo, stla, stlo, flattening, &azimuth,
                                &back_azimuth);
    if (dist<0.001) // When dist is small, using sac to read the file will result in N/A dist.
        out_CC.shd.lcalda=0;
    out_CC.shd.gcarc=dist*180/PI;
    out_CC.shd.dist=dist*semimajor;
    out_CC.shd.az= azimuth*180./PI;
    out_CC.shd.baz= back_azimuth*180./PI-180.0;
    out_CC.shd.user0=1.0; // Stack day
    out_CC.shd.user1=sac_am1.shd.cmpaz; // for component correction
    out_CC.shd.user2=sac_am2.shd.cmpaz;
    /*--- end of setting header value ----*/
    return true;
}


CC_output::CC_output(const std::string & OUTname, const std::string & MDIR, const std::string  & STAN1, const std::string  & CHAN1,
                     const std::string  & STAN2, const std::string  & CHAN2, const int & DAY)
{
    outfname.push_back(OUTname);
    day.push_back(DAY);
    monthdir=MDIR;
    sta1=STAN1;
    chan1=CHAN1;
    sta2=STAN2;
    chan2=CHAN2;
}
CC_output::CC_output( ) {}
CC_output::~CC_output( ) {}

CC_pairs::CC_pairs() {}
CC_pairs::~CC_pairs() {}
CC_pairs::CC_pairs(const std::string  & STAN1, const std::string  & CHAN1,
                   const std::string  & STAN2, const std::string  & CHAN2)
{
    sta1=STAN1;
    sta2=STAN2;
    chan1=CHAN1;
    chan2=CHAN2;
}

sta_pairs::sta_pairs() {}
sta_pairs::~sta_pairs() {}
sta_pairs::sta_pairs(const std::string  & STAN1, const std::string  & STAN2)
{
    sta1=STAN1;
    sta2=STAN2;
}

Station::Station() {}
Station::~Station() {}
Station::Station(const std::string  & STAN,  const int & CCFLAG)
{
    sta=STAN;
    CCflag=CCFLAG;
}

bool Station::checkdoCC(const Station & STA )
{
    if  (CCflag==STA.CCflag && CCflag==0)
        return false;
    else if(CCflag!=STA.CCflag && CCflag!=0 && STA.CCflag!=0)
        return false;
    else if (CCflag< 0 && sta==STA.sta)
    {
        //logger.Hold( INFO, "GROUP: "+ std::to_string(groupflag)+". "+"DO NOT DO AUTO", FuncName );
        std::cout << "GROUP: "<<CCflag<<std::endl;
        std::cout<<"DO NOT DO AUTO"<<std::endl;
        return false;
    }
    else
        return true;
}

outCCRec::outCCRec() {}
outCCRec::~outCCRec() {}
outCCRec::outCCRec(std::vector < SacRec > & INRec, int CCSIZE, int DAY, int LISTNo, std::string MDIR,
                   std::string STA1, std::string STA2, std::vector <std::string> & OUTFNAME, std::string CORDIR)
{
    for (int i=0; i<CCSIZE; i++)
    {
        if ( INRec[i].sig ) CCRec.push_back( std::move(INRec[i]) );
        outfname.push_back(std::move(OUTFNAME[i]));
    }
    sacRecSize=CCSIZE;
    day=DAY;
    cclistNo=LISTNo;
    monthdir=MDIR;
    sta1=STA1;
    sta2=STA2;
    COR=CORDIR;
}

bool Rotation(sta_pairs & stapair, std::vector < std::string > & CHAN, std::vector < std::string > & DIR, std::string CHPRE)
{
    if (CHAN.size()<2)
    {
        logger.Hold( WARNING, "Only One component in channel list, cannot do rotation!", FuncName );
        //std::cout<<"Waring: Only One component in channel list, cannot do rotation!"<<std::endl;
        return false;
    }
    //int rotflag=0;
    float psi, theta,Cpsi,Spsi,Ctheta,Stheta;
    const float PI=3.14159265358979323846;
    int N;
    for (int i=0; i<DIR.size(); i++)
    {
        std::string fEE,fEN,fNN,fNE,fEZ,fZE,fNZ,fZN;
        std::string sta1=stapair.sta1;
        std::string sta2=stapair.sta2;
        std::string dir=DIR[i];
        if (sta1>sta2 || sta1== sta2) continue;
        fEE=dir+"/"+sta1+"/COR_"+sta1+"_"+CHAN[0]+"_"+sta2+"_"+CHAN[0]+".SAC";
        fNN=dir+"/"+sta1+"/COR_"+sta1+"_"+CHAN[1]+"_"+sta2+"_"+CHAN[1]+".SAC";
        fEN=dir+"/"+sta1+"/COR_"+sta1+"_"+CHAN[0]+"_"+sta2+"_"+CHAN[1]+".SAC";
        fNE=dir+"/"+sta1+"/COR_"+sta1+"_"+CHAN[1]+"_"+sta2+"_"+CHAN[0]+".SAC";

        if ( FileExists(fEE)&&FileExists(fNN)&&FileExists(fEN)&&FileExists(fNE) )
        {
            //std::cout<<"Do Rotation(RT) for:"<<sta1<<" and "<<sta2<<" at dir:"<<dir<<std::endl;
            logger.Hold( INFO, "Do Rotation(RT) for:"+sta1+" and "+sta2+" at dir:"+dir, FuncName );
            SacRec CorEE,CorEN,CorNN,CorNE;
            SacRec CorTT,CorRR,CorTR,CorRT;
            CorEE.Load(fEE);
            CorEN.Load(fEN);
            CorNN.Load(fNN);
            CorNE.Load(fNE);
            CorTT.Load(fEN);
            CorRR.Load(fEN);
            CorTR.Load(fEN);
            CorRT.Load(fEN);
            theta=CorEN.shd.az;
            psi=CorEN.shd.baz;
            //  std::cout<<"THETA: "<<theta<<" PSI: "<<psi<<std::endl;
            Ctheta=cos(PI*theta/180);
            Stheta=sin(PI*theta/180);
            Cpsi=cos(PI*psi/180);
            Spsi=sin(PI*psi/180);
            // std::cout<<"COS THETA: "<<Ctheta<<" COS PSI: "<<Cpsi<<std::endl;
            // std::cout<<"SIN THETA: "<<Stheta<<" SIN PSI: "<<Spsi<<std::endl;
            N=CorEE.shd.npts;
            for (int r =0; r < N; r++)
            {
                CorTT.sig[r] = -Ctheta*Cpsi*CorEE.sig[r] +
                               Ctheta*Spsi*CorEN.sig[r] - Stheta*Spsi*CorNN.sig[r] + Stheta*Cpsi*CorNE.sig[r];

                CorRR.sig[r] = - Stheta*Spsi*CorEE.sig[r] -
                               Stheta*Cpsi*CorEN.sig[r] - Ctheta*Cpsi*CorNN.sig[r] - Ctheta*Spsi*CorNE.sig[r];

                CorTR.sig[r] = -Ctheta*Spsi*CorEE.sig[r] -
                               Ctheta*Cpsi*CorEN.sig[r] + Stheta*Cpsi*CorNN.sig[r] + Stheta*Spsi*CorNE.sig[r];

                CorRT.sig[r] = -Stheta*Cpsi*CorEE.sig[r] +
                               Stheta*Spsi*CorEN.sig[r] + Ctheta*Spsi*CorNN.sig[r] - Ctheta*Cpsi*CorNE.sig[r];
            }

            strcpy(CorTT.shd.kcmpnm,(CHPRE+"T"+CHPRE+"T").c_str());
            strcpy(CorRR.shd.kcmpnm,(CHPRE+"R"+CHPRE+"R").c_str());
            strcpy(CorTR.shd.kcmpnm,(CHPRE+"T"+CHPRE+"R").c_str());
            strcpy(CorRT.shd.kcmpnm,(CHPRE+"R"+CHPRE+"T").c_str());
            CorTT.Write(dir+"/"+sta1+"/COR_"+sta1+"_"+CHPRE+"T"+"_"+sta2+"_"+CHPRE+"T"+".SAC");
            CorRR.Write(dir+"/"+sta1+"/COR_"+sta1+"_"+CHPRE+"R"+"_"+sta2+"_"+CHPRE+"R"+".SAC");
            CorTR.Write(dir+"/"+sta1+"/COR_"+sta1+"_"+CHPRE+"T"+"_"+sta2+"_"+CHPRE+"R"+".SAC");
            CorRT.Write(dir+"/"+sta1+"/COR_"+sta1+"_"+CHPRE+"R"+"_"+sta2+"_"+CHPRE+"T"+".SAC");
        }
        if (CHAN.size()==2)
            continue;
        // Rotation of RTZ
        fEZ=dir+"/"+sta1+"/COR_"+sta1+"_"+CHAN[0]+"_"+sta2+"_"+CHAN[2]+".SAC";
        fZE=dir+"/"+sta1+"/COR_"+sta1+"_"+CHAN[2]+"_"+sta2+"_"+CHAN[0]+".SAC";
        fNZ=dir+"/"+sta1+"/COR_"+sta1+"_"+CHAN[1]+"_"+sta2+"_"+CHAN[2]+".SAC";
        fZN=dir+"/"+sta1+"/COR_"+sta1+"_"+CHAN[2]+"_"+sta2+"_"+CHAN[1]+".SAC";
        if ( FileExists(fEZ)&&FileExists(fZE)&&FileExists(fNZ)&&FileExists(fZN) )
        {
            logger.Hold( INFO, "Do Rotation(RTZ) for:"+sta1+" and "+sta2+" at dir:"+dir, FuncName );
            //std::cout<<"Do Rotation(RTZ) for:"<<sta1<<" and "<<sta2<<" at dir:"<<dir<<std::endl;
            SacRec CorEZ,CorZE,CorNZ,CorZN;
            SacRec CorRZ,CorZR,CorTZ,CorZT;
            CorEZ.Load(fEZ);
            CorZE.Load(fZE);
            CorNZ.Load(fNZ);
            CorZN.Load(fZN);
            CorRZ.Load(fEZ);
            CorZR.Load(fEZ);
            CorTZ.Load(fEZ);
            CorZT.Load(fEZ);
            theta=CorEZ.shd.az;
            psi=CorEZ.shd.baz;
            //  std::cout<<"THETA: "<<theta<<" PSI: "<<psi<<std::endl;
            Ctheta=cos(PI*theta/180);
            Stheta=sin(PI*theta/180);
            Cpsi=cos(PI*psi/180);
            Spsi=sin(PI*psi/180);
            // std::cout<<"COS THETA: "<<Ctheta<<" COS PSI: "<<Cpsi<<std::endl;
            //  std::cout<<"SIN THETA: "<<Stheta<<" SIN PSI: "<<Spsi<<std::endl;
            N=CorEZ.shd.npts;
            for ( int r =0; r < N; r++ )
            {
                CorRZ.sig[r] = Ctheta*CorNZ.sig[r] + Stheta*CorEZ.sig[r];

                CorZR.sig[r] = - Cpsi*CorZN.sig[r] -Spsi*CorZE.sig[r];

                CorTZ.sig[r] = -Stheta*CorNZ.sig[r] + Ctheta*CorEZ.sig[r];

                CorZT.sig[r] =  Spsi*CorZN.sig[r] - Cpsi*CorZE.sig[r];
            }
            strcpy(CorRZ.shd.kcmpnm,(CHPRE+"R"+CHPRE+"Z").c_str());
            strcpy(CorZR.shd.kcmpnm,(CHPRE+"Z"+CHPRE+"R").c_str());
            strcpy(CorTZ.shd.kcmpnm,(CHPRE+"T"+CHPRE+"Z").c_str());
            strcpy(CorZT.shd.kcmpnm,(CHPRE+"Z"+CHPRE+"T").c_str());
            CorRZ.Write(dir+"/"+sta1+"/COR_"+sta1+"_"+CHPRE+"R"+"_"+sta2+"_"+CHPRE+"Z"+".SAC");
            CorZR.Write(dir+"/"+sta1+"/COR_"+sta1+"_"+CHPRE+"Z"+"_"+sta2+"_"+CHPRE+"R"+".SAC");
            CorTZ.Write(dir+"/"+sta1+"/COR_"+sta1+"_"+CHPRE+"T"+"_"+sta2+"_"+CHPRE+"Z"+".SAC");
            CorZT.Write(dir+"/"+sta1+"/COR_"+sta1+"_"+CHPRE+"Z"+"_"+sta2+"_"+CHPRE+"T"+".SAC");
        } //End for loop sta2
    } // End of dir
    return true;
}

void StackAll(std::vector < Station > & StaList, std::vector < std::string > & MDIR, std::vector < std::string > & CHAN, CCDatabase & CDB)
{
    if (CDB.GetParams().fstackall!=0)
    {
        std::cout<<"ATTENTION: Start Stacking All "<<std::endl;
        std::vector  <CC_pairs> pair_list;
        for (int s1=0; s1<StaList.size(); s1++)
            for (int s2=0; s2<StaList.size(); s2++)
            {
                if (StaList[s1].sta>StaList[s2].sta || !StaList[s1].checkdoCC(StaList[s2])) continue;
                for (int ch1=0; ch1<CHAN.size(); ch1++)
                    for (int ch2=0; ch2<CHAN.size(); ch2++)
                        pair_list.push_back(CC_pairs(StaList[s1].sta, CHAN[ch1], StaList[s2].sta, CHAN[ch2]));
            }
        bool delete_monflag=false;
        if (CDB.GetParams().fstackall==2) delete_monflag=true;
        #pragma omp parallel for
        for (int i=0; i<pair_list.size(); i++)
        {
            CC_output all_stack;
            bool newflag=true;
            for (int j = 0; j < MDIR.size(); j++)
            {
                std::string monthly_name=MDIR[j]+"/COR/"+pair_list[i].sta1+"/COR_"+pair_list[i].sta1+
                                         "_"+pair_list[i].chan1+"_"+pair_list[i].sta2+"_"+pair_list[i].chan2+".SAC";
                if (!FileExists(monthly_name)) continue;
                if(newflag==true)
                {
                    all_stack=CC_output(monthly_name, MDIR[j], pair_list[i].sta1, pair_list[i].chan1, pair_list[i].sta2, pair_list[i].chan2, 0);
                    newflag=false;
                    continue;
                }
                all_stack.outfname.push_back(monthly_name);
                all_stack.day.push_back(0);
            }
            if (all_stack.day.empty()) continue;
            stackCCdata(all_stack,  0);
            if (delete_monflag==true)
                dRemove((all_stack.monthdir+"/COR").c_str());
        }
        std::cout<<"ATTENTION: End of  Stacking All "<<std::endl;
    }
}


bool CalcRecCor( std::string & fname1, std::string fname2, float *cor_rec, int lagn, int mintlen, float tlen, float fs)
{
    int t, irec1, irec2;
    int recB, recE;
    int TLEN=(int)(tlen*fs)+1;
    int recbeg, recend,n;
    FILE *fin1, *fin2;
    std::vector <hole_rec> sta1, sta2;
    //std::cout<<"f1 rec: "<<fname1<<std::endl;
    if (FileExists(fname1))
    {
        if ( (fin1=fopen((fname1).c_str(),"r")) !=NULL );
        while((n = fscanf(fin1,"%d %d", &recbeg, &recend)) != EOF) // May cause error if type does not match exactly!!!
        {
            sta1.push_back(hole_rec(recbeg,recend));
        }
    }
    else
        sta1.push_back(hole_rec(0, TLEN-1));
    if (FileExists(fname2))
    {
        if ( (fin2=fopen((fname2).c_str(),"r")) !=NULL );
        while((n = fscanf(fin2,"%d %d", &recbeg, &recend)) != EOF) // May cause error if type does not match exactly!!!
        {
            sta2.push_back(hole_rec(recbeg,recend));
        }
    }
    else
        sta2.push_back(hole_rec(0, TLEN-1));
    for(t=0; t<=lagn; t++)
    {
        cor_rec[lagn+t]=1;
        cor_rec[lagn-t]=1;
        for(irec1=0; irec1<sta1.size(); irec1++)
        {
            for(irec2=0; irec2<sta2.size(); irec2++)
            {
                if(sta1[irec1].rec_b >=sta2[irec2].rec_e-t) continue;
                if(sta1[irec1].rec_e<=sta2[irec2].rec_b-t) break;
                recB = std::max(sta1[irec1].rec_b, sta2[irec2].rec_b-t);
                recE = std::min(sta1[irec1].rec_e, sta2[irec2].rec_e-t);
                cor_rec[lagn+t] += recE - recB;
            }
            for(irec2=0; irec2<sta2.size(); irec2++)
            {
                if(sta1[irec1].rec_b >=sta2[irec2].rec_e+t) continue;
                if(sta1[irec1].rec_e<=sta2[irec2].rec_b+t) break;
                recB = std::max(sta1[irec1].rec_b, sta2[irec2].rec_b+t);
                recE = std::min(sta1[irec1].rec_e, sta2[irec2].rec_e+t);
                cor_rec[lagn-t] += recE - recB;
            }
        }
    }
    cor_rec[lagn] /= 2;

    if(cor_rec[0]< (int)(mintlen*fs) || cor_rec[lagn*2]< (int)(mintlen*fs) )
    {
        std::cout<<"*** cor time less than "<<mintlen<< "sec. Skipped! *** " <<std::endl;
        return false;
    }
    return true;
}

hole_rec::hole_rec() {}
hole_rec::~hole_rec() {}
hole_rec::hole_rec(int BEG, int END)
{
    rec_b=BEG;
    rec_e=END;

}
