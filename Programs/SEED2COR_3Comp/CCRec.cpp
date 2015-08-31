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

bool FileExists(const char* filename)
{
    struct stat info;
    int ret = -1;

    //get the file attributes
    ret = stat(filename, &info);
    if(ret == 0)
    {
        return true;
    }
    else
    {
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


bool doCor(SacRec & out_CC, const SacRec & sac_am1, const SacRec & sac_ph1, const SacRec & sac_am2, const SacRec & sac_ph2,
           const int & lagtime, const int & fs, const int & mintlen, const float & tlen, const int & ftlen, std::string & fname1, std::string & fname2,
           const double dist, const double azimuth, const double back_azimuth)  // Note that lagtime and fs are integer!!!
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
    const double semimajor = 6378.137;
    const double PI = 3.14159265358979323846;
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

outCCRec::outCCRec() {}
outCCRec::~outCCRec() {}
outCCRec::outCCRec(std::vector < SacRec > & INRec, int CCSIZE, int DAY, int LISTNo, std::string MDIR,
             std::string STA1,  std::vector < std::string >  ChanL1, std::string STA2, std::vector < std::string >  ChanL2,
             std::vector <std::string> & OUTFNAME, std::string CORDIR)
{
    for (int i=0; i<CCSIZE; i++)
    {
        if ( INRec[i].sig ) CCRec.push_back( std::move(INRec[i]) );
        outfname.push_back(std::move(OUTFNAME[i]));
    }
        for (int i=0; i<ChanL1.size(); i++)
    {
        chanL1.push_back(std::move(ChanL1[i]));
    }
        for (int i=0; i<ChanL2.size(); i++)
    {
        chanL2.push_back(std::move(ChanL2[i]));
    }
    sacRecSize=CCSIZE;
    day=DAY;
    cclistNo=LISTNo;
    monthdir=MDIR;
    sta1=STA1;
    sta2=STA2;
    COR=CORDIR;
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
