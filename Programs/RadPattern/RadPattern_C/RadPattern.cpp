#include "RadPattern.h"
#include <fstream>
#include <cstring>
#include <algorithm>

/* FORTRAN entrance */
const int nazi = RadPattern::nazi;
extern"C" {
   void rad_pattern_r_(char *feig_buff, int *eig_len, int *phvnper, float *phvdper,
                       const float *strike, const float *dip, const float *rake, const float *depth,
							  const float *per, int *nper,
                       float *azi, float grT[][nazi], float phT[][nazi], float amp[][nazi]);
   void rad_pattern_l_(char *feig_buff, int *eig_len, int *phvnper, float *phvdper,
                       const float *strike, const float *dip, const float *rake, const float *depth,
							  const float *per, int *nper,
                       float *azi, float grT[][nazi], float phT[][nazi], float amp[][nazi]);
}


/* ---------- implimentation ---------- */
struct RadPattern::Rimpl {
   //FocalInfo<ftype> finfo, finfoold;
   //std::string outname_mis;
   std::string feignmem, fphvnmem; // name of files of the currently-in-the-memory contents

	static const int NaN = RadPattern::NaN;
   int feig_len = NaN;
   int phvnper = NaN; float phvdper = NaN;
   char *feig_buff = nullptr;

   /* ---------- con/destructors ---------- */
	Rimpl() {}

	Rimpl( const Rimpl& r2 ) 
		: phvnper(r2.phvnper), phvdper(r2.phvdper)
		, feig_len(r2.feig_len) {
		if( feig_len > 0) {
			feig_buff = new char[feig_len];
			memcpy(feig_buff, r2.feig_buff, feig_len );
		}
	}

   ~Rimpl() {
		if(feig_buff) {
			delete [] feig_buff;
			feig_buff = nullptr;
		}
	}

};


/* con/destructors and operators */
RadPattern::RadPattern()
   : pimplR( new Rimpl ) {}

RadPattern::RadPattern( const RadPattern& rp2 )
   : pimplR( new Rimpl(*(rp2.pimplR)) ), type( rp2.type )
   , stk( rp2.stk ), dip( rp2.dip )
   , rak( rp2.rak ), dep( rp2.dep )
   , aziV( rp2.aziV ), grtM( rp2.grtM ) 
   , phtM( rp2.phtM ), ampM( rp2.ampM ) {}

RadPattern::RadPattern( RadPattern&& rp2 )
   : pimplR( std::move(rp2.pimplR) ), type( rp2.type )
	, stk( rp2.stk ), dip( rp2.dip )
	, rak( rp2.rak ), dep( rp2.dep )
	, aziV( std::move(rp2.aziV) ), grtM( std::move(rp2.grtM) )
	, phtM( std::move(rp2.phtM) ), ampM( std::move(rp2.ampM) ) {}

RadPattern& RadPattern::operator= ( const RadPattern& rp2 ) {
   pimplR.reset( new Rimpl(*(rp2.pimplR)) );
	type = rp2.type;
	stk = rp2.stk; dip = rp2.dip;
	rak = rp2.rak; dep = rp2.dep;
	aziV = rp2.aziV; grtM = rp2.grtM;
	phtM = rp2.phtM; ampM = rp2.ampM;
}

RadPattern& RadPattern::operator= ( RadPattern&& rp2 ){
   pimplR = std::move(rp2.pimplR);
	type = rp2.type;
	stk = rp2.stk; dip = rp2.dip;
	rak = rp2.rak; dep = rp2.dep;
	aziV = std::move(rp2.aziV); grtM = std::move(rp2.grtM);
	phtM = std::move(rp2.phtM); ampM = std::move(rp2.ampM);
}

RadPattern::~RadPattern() {}

/* copy arrayin[nazi] into Vout with positions shifted by int(nazi/2) */
void RadPattern::ShiftCopy( std::vector<float>& Vout, const float* arrayin, const int nazi ) const {
   // check nazi
   if( nazi % 2 == 0 )
      throw ErrorRP::BadParam( FuncName, "unexpected nazi = " + std::to_string(nazi) );
   int nazio2 = nazi / 2;
   Vout.clear(); Vout.reserve(nazi);
   // shift by 180 degree
   Vout = std::vector<float>( arrayin+nazio2, arrayin+nazi-1 );
   Vout.insert( Vout.end(), arrayin, arrayin+nazio2+1 );
}

/* predict radpattern for rayleigh and love waves */
bool RadPattern::Predict( char typein, const std::string& feigname, const std::string& fphvname,
			  const ftype stkin, const ftype dipin, const ftype rakin, const ftype depin,
			  const std::vector<float>& perlst ) {
			  //std::vector< std::vector<AziData> >& per_azi_pred ) {
   if( typein!='R' && typein!='L' )
      throw ErrorRP::BadParam(FuncName, "unknown type = "+type);

	// return if the requested new state is exactly the same as the one stored
	if( type==typein && stk==stkin && dip==dipin && 
		 rak==rakin && dep==depin && perlst.size()<=grtM.size() ) {
		bool allfound = true;
		for( const auto per : perlst )
			if( grtM.find(per) == grtM.end() ) {
				allfound = false;
				break;
			}
		if( allfound ) return false;	// not updated
	}

	// store current state;
	type = typein; 
	stk = stkin; dip = dipin;
	rak = rakin; dep = depin;

  #pragma omp critical
  { // omp critical begins
   // read feig into memory
   if( pimplR->feig_buff == nullptr || feigname != pimplR->feignmem ) {
      if( feigname.empty() ) throw ErrorRP::BadParam(FuncName, "empty feigname");
      std::ifstream fin( feigname.c_str() );
      if( ! fin ) throw ErrorRP::BadFile(FuncName, feigname);
      fin.seekg(0, std::ios::end);
      pimplR->feig_len = fin.tellg();
      if( pimplR->feig_buff ) {
			delete [] pimplR->feig_buff;
			pimplR->feig_buff = nullptr;
      }
      pimplR->feig_buff = new char[pimplR->feig_len];
      fin.seekg(0,std::ios::beg);
      fin.read(pimplR->feig_buff, pimplR->feig_len);
      fin.close();
      pimplR->feignmem = feigname;
   }

   // read in nper and dper from fphv
   if( pimplR->phvnper < 0 || fphvname != pimplR->fphvnmem )  {
      if( fphvname.empty() ) throw ErrorRP::BadParam(FuncName, "empty fphvname");
      std::ifstream fin( fphvname.c_str() );
      if( ! fin ) throw ErrorRP::BadFile(FuncName, fphvname);
      pimplR->phvnper = 0;
      for( std::string line; std::getline(fin, line); ) {
         float ftmp1, ftmp2, ftmp3;
         if( sscanf(line.c_str(), "%f %f %f", &ftmp1, &ftmp2, &ftmp3) != 3 ) continue;
         if( pimplR->phvnper == 0 ) pimplR->phvdper = ftmp1; //std::cerr<<"dper = "<<pimplR->phvdper<<" ftmp1 = "<<ftmp1<<std::endl; }
         else if( pimplR->phvnper == 1 ) pimplR->phvdper = ftmp1 - pimplR->phvdper; //std::cerr<<"dper = "<<pimplR->phvdper<<" ftmp1 = "<<ftmp1<<std::endl; }
         (pimplR->phvnper)++;
      }
      fin.close();
      pimplR->fphvnmem = fphvname;
   }
  } // omp critical ends

   // check if feig and fphv contents are modified
   if( pimplR->feig_buff == nullptr || feigname != pimplR->feignmem ) 
      throw ErrorRP::BadBuff(FuncName, feigname + " != " + pimplR->feignmem);
   if( pimplR->phvnper < 0 || fphvname != pimplR->fphvnmem ) 
      throw ErrorRP::BadBuff(FuncName, fphvname + " != " + pimplR->fphvnmem);

   // run rad_pattern
   int nper = perlst.size();
	float azi[nazi], grT[nper][nazi], phT[nper][nazi], amp[nper][nazi];

   if( type == 'R' ) {
      rad_pattern_r_( pimplR->feig_buff, &(pimplR->feig_len), &(pimplR->phvnper), &(pimplR->phvdper),
                      &(stk), &(dip), &(rak), &(dep), &(perlst.at(0)), &nper, azi, grT, phT, amp );
      //std::cerr<<pimplR->feignmem<<" "<<pimplR->feig_len<<" "<<pimplR->phvnper<<" "<<pimplR->phvdper<<"  "<<stk<<" "<<dip<<" "<<rak<<" "<<dep<<"   "<<nper<<std::endl;
   } else if( type == 'L' ) {
      rad_pattern_l_( pimplR->feig_buff, &(pimplR->feig_len), &(pimplR->phvnper), &(pimplR->phvdper),
                      &(stk), &(dip), &(rak), &(dep), &(perlst.at(0)), &nper, azi, grT, phT, amp );
   }

   // check if feig and fphv contents are modified
   if( pimplR->feig_buff == nullptr || feigname != pimplR->feignmem ) 
      throw ErrorRP::BadBuff(FuncName, feigname + " != " + pimplR->feignmem);
   if( pimplR->phvnper < 0 || fphvname != pimplR->fphvnmem ) 
      throw ErrorRP::BadBuff(FuncName, fphvname + " != " + pimplR->fphvnmem);

	// copy predictions into maps
	grtM.clear(); phtM.clear(); ampM.clear(); //aziV.clear();
   //float azi[nazi], grT[nper][nazi], phT[nper][nazi], amp[nper][nazi];
	for( int iper=0; iper<perlst.size(); iper++ ) {
		float per = perlst[iper];

/*
		// copy as is
		auto &grV = grtM[per], &phV = phtM[per], &amV = ampM[per];
      grV = std::vector<float>( grT[iper], grT[iper]+nazi );
      phV = std::vector<float>( phT[iper], phT[iper]+nazi );
      amV = std::vector<float>( amp[iper], amp[iper]+nazi );
		// shift the phase by pi
		float T_pi = per * 0.5;
		auto pi_shifter = [&](float& val) {
			val += T_pi;
			if( val > T_pi ) val -= per;
		};
		std::for_each( phV.begin(), phV.end(), pi_shifter );
		// flip
		std::transform( grV.begin(), grV.end(), grV.begin(), std::negate<float>() );
		std::transform( phV.begin(), phV.end(), phV.begin(), std::negate<float>() );
*/
		// shift by 180 degree
		ShiftCopy( grtM[per], grT[iper], nazi );
		ShiftCopy( phtM[per], phT[iper], nazi );
		ShiftCopy( ampM[per], amp[iper], nazi );

   }
	//aziV = std::vector<float>( azi, azi+nazi );

	// invalidate focal predictions with amplitudes < amp_avg * AmpValidPerc
	//const float NaN = Rimpl::NaN;
	for( int iper=0; iper<perlst.size(); iper++ ) {
		float per = perlst[iper];
		// determine min amplitude
		int Nvalid = 0; float Amin = 0.;
		const auto& ampV = ampM[per];
		for( const auto& amp : ampV ) Amin += amp;
		Amin *= (AmpValidPerc / ampV.size());
		// invalidate azimuths with small amplitudes
		auto& grtV = grtM[per];
		for(int iazi=0; iazi<nazi; iazi++)
			if( ampV[iazi] < Amin ) {
				int jazilow = iazi-InvalidateHwidth, jazihigh = iazi+InvalidateHwidth+1;
				// invalidate the range (jazilow - 0)
				int jazi = jazilow;
				if( jazi < 0 )
					for(; jazi<0; jazi++)
						grtV[jazi + nazi] = NaN;
				// invalidate the range (0 - 360)
				for(; jazi<jazihigh&&jazi<nazi; jazi++)
					grtV[jazi] = NaN;
				// invalidate the range (360 - jazihigh)
				if( jazihigh > nazi )
					for(; jazi<jazihigh; jazi++)
						grtV[jazi - nazi] = NaN;
			}
	}

	return true;	// updated!
}


bool RadPattern::GetPred( const float per, const float azi,
								  float& grt, float& pht, float& amp ) const {
	// check validities of period and azimuth
	auto Igrt = grtM.find(per);
	if( Igrt == grtM.end() )
		throw ErrorRP::BadParam(FuncName, "un-predicted period");
	int iazi = (int)(azi/dazi);
	if( iazi<0 || iazi >= nazi )
		throw ErrorRP::BadAzi( FuncName, "azi = "+std::to_string(azi) );
	// low and high azimuth
	float azil = iazi*dazi, azih = azil+dazi;
	float azifactor = (azi-azil) / dazi, ftmp1, ftmp2;
	//const float NaN = Rimpl::NaN;
	// group delay
	ftmp1 = (Igrt->second)[iazi], ftmp2 = (Igrt->second)[iazi+1];
	if( ftmp1==NaN || ftmp2==NaN ) {
		grt = pht = amp = NaN;
		return false;
	}
	grt = ftmp1 + (ftmp2 - ftmp1) * azifactor;
	// phase shift
	auto Ipht = phtM.find(per);
	ftmp1 = (Ipht->second)[iazi], ftmp2 = (Ipht->second)[iazi+1];
	pht = ftmp1 + (ftmp2 - ftmp1) * azifactor;
	// amp
	auto Iamp = ampM.find(per);
	ftmp1 = (Iamp->second)[iazi], ftmp2 = (Iamp->second)[iazi+1];
	amp = ftmp1 + (ftmp2 - ftmp1) * azifactor;

	return true;
}

void RadPattern::OutputPreds( const std::string& fname, const float Afactor ) {
	if( grtM.size() == 0 ) return;

	std::ofstream fout( fname );
   if( ! fout ) throw ErrorRP::BadFile(FuncName, fname);
	for( auto grtP : grtM ) {
		float per = grtP.first;
		auto grtV = grtP.second;
		auto phtV = phtM.find(per)->second;
		auto ampV = ampM.find(per)->second;
		for( int iazi=0; iazi<nazi; iazi++ ) {
			float azi = iazi*dazi;
			float grt = grtV.at(iazi);
			if( grt == RadPattern::NaN ) continue;
			float pht = phtV.at(iazi);
			float amp = ampV.at(iazi);
			fout<<azi<<" "<<grt<<" "<<pht<<" "<<amp*Afactor<<" "<<per<<"\n";
			//std::cerr<<azi<<" "<<grt<<" "<<pht<<" "<<amp*Afactor<<" "<<per<<"\n";
		}
		fout<<"\n\n";
	}
}
