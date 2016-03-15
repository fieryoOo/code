#include "SacRec.h"
#include "InfoLists.h"
#include <iostream>
#include <vector>
#include <deque>
#include <cmath>

bool FileExists(const char* filename);
bool FileExists(const std::string& filename);

/* normalize all sac files in sacV (simultaneously if SyncNorm==true) */
void TNormAll( std::deque<SacRec>& sacV, const std::vector<DailyInfo>& dinfoV, bool SyncNorm)
{
    if( sacV.empty() ) return;
    if( sacV.size()!=dinfoV.size() )
        throw std::runtime_error("Error(TNormAll): size mismatch ("+std::to_string(sacV.size())+" - "+std::to_string(dinfoV.size()));
    SacRec sac_sigmax;
    for( int isac=0; isac<sacV.size(); isac++ )
    {
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
                    //sac.Filter( f1, f2, f3, f4, sac_eqk );
						  sac.BandpassCOSFilt( f1, f2, f3, f4, sac_eqk );
                    sac_eqk.Smooth( dinfo.timehlen, sac_sm );
                } else {
                    sac.Smooth( dinfo.timehlen, sac_sm );
                }
					// synchronize sac_sm
					float sacmean, sacstd;
					sac_sm.MeanStd ( -12345.,  -12345., sacmean, sacstd );
					sacmean = 1./fabs(sacmean);
					sac_sm.Mul(sacmean);
					sac.shd.user1 = sacmean;
                /* max smoothed signal */
                sac_sigmax.PullUpTo( sac_sm );
            } else {
                sac.RunAvg( dinfo.timehlen, dinfo.Eperl, dinfo.Eperh );
            }
            break;
        case 3:
            throw std::runtime_error("EqkCut not exist!");/// TODO: add earthquake cutting normalization
            break;
        }
    }
    if( SyncNorm  && dinfoV[0].tnorm_flag==2 ) /// 1.02-->1.03
        for( auto& sac : sacV )
            if( sac.sig ) {
					sac.Divf( sac_sigmax );
					sac.Mul(sac.shd.user1);
				}
    //sacV[0].Write( "Test_RunAvg.sac" );
}

/* normalize and apply taper */
void FNormAll( std::deque<SacRec>& sacV, const std::vector<DailyInfo>& dinfoV, bool SyncNorm )
{
    if( sacV.empty() ) return;
    if( sacV.size()!=dinfoV.size() )
        throw std::runtime_error("Error(FNormAll): size mismatch ("+std::to_string(sacV.size())+" - "+std::to_string(dinfoV.size()));

    // compute max smoothed signal
    SacRec sac_sigmax;
    for( int isac=0; isac<sacV.size(); isac++ )
    {
        auto& dinfo = dinfoV[isac];
        auto& sac = sacV[isac];
        if( !(sac.sig) ) continue;
        if( SyncNorm )
        {
            SacRec sac_sm;
            sac.Smooth( dinfo.frechlen, sac_sm );
            sac_sigmax.PullUpTo( sac_sm );
        }
        else
        {
            sac.RunAvg( dinfo.frechlen, -1., -1. );
        }
    }
    // normalize and taper
    if( SyncNorm )
        for( int isac=0; isac<sacV.size(); isac++ )
        {
            auto& dinfo = dinfoV[isac];
            auto& sac = sacV[isac];
            if( ! sac.sig ) continue;
            sac.Divf( sac_sigmax );
            float fl=1./dinfo.perh, fh=1./dinfo.perl;
            sac.cosTaperL( fl*0.8, fl );
            sac.cosTaperR( fh, fh*1.2 );
        }
}
