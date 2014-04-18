#include "RadPattern.h"
#include <fstream>

/* FORTRAN entrance */
extern"C" {
   void rad_pattern_r_(char *feig_buff, int *eig_len, int *phvnper, int *phvdper,
                       float *strike, float *dip, float *rake, float *depth, float *per, int *nper,
                       float *azi, float grT[][181], float phT[][181], float amp[][181]);
   void rad_pattern_l_(char *feig_buff, int *eig_len, int *phvnper, int *phvdper,
                       float *strike, float *dip, float *rake, float *depth, float *per, int *nper,
                       float *azi, float grT[][181], float phT[][181], float amp[][181]);
}


/* implimentation */
struct RadPattern::Rimpl {
   FocalInfo<float> finfo, finfoold;
   std::string outname_mis;
   std::string feignmem, fphvnmem; // name of files of the currently-in-the-memory contents

   int phvnper, phvdper;
   int feig_len;
   char *feig_buff;

   Rimpl() {
      phvnper = -12345; phvdper = -12345;
      feig_buff = NULL; feig_len = -12345;
   }
   ~Rimpl() { if(feig_buff) free(feig_buff); }

};


/* con/destructors */
RadPattern::RadPattern()
   : pimplR(new Rimpl) {}

RadPattern::~RadPattern() {}

/* predict radpattern for rayleigh and love waves */
bool RadPattern::Predict( char type, const std::string& feigname, const std::string& fphvname, const FocalInfo<float>& finfo,
				 std::vector<float>& perlst, std::vector< std::vector<AziData> >& per_azi_pred ) {
   if( type!='R' && type!='L' ) return false;
   int nper = perlst.size();
   float strike = finfo.strike, dip = finfo.dip, rake = finfo.rake, depth = finfo.depth;
   float azi[181], grT[nper][181], phT[nper][181], amp[nper][181];

   // read feig into memory
   if( pimplR->feig_buff == NULL || feigname != pimplR->feignmem ) {
      if( feigname.empty() ) return false;
      std::ifstream fin( feigname.c_str() );
      if( ! fin ) return false;
      fin.seekg(0, std::ios::end);
      pimplR->feig_len = fin.tellg();
      pimplR->feig_buff = new char[pimplR->feig_len];
      fin.seekg(0,std::ios::beg);
      fin.read(pimplR->feig_buff, pimplR->feig_len);
      fin.close();
      pimplR->feignmem = feigname;
   }

   // read in nper and dper from fphv
   if( pimplR->phvnper < 0 || fphvname != pimplR->fphvnmem )  {
      if( fphvname.empty() ) return false;
      std::ifstream fin( fphvname.c_str() );
      if( ! fin ) return false;
      pimplR->phvnper = 0;
      for( std::string line; std::getline(fin, line); ) {
         float ftmp1, ftmp2, ftmp3;
         if( sscanf(line.c_str(), "%f %f %f", &ftmp1, &ftmp2, &ftmp3) != 3 ) continue;
         if( pimplR->phvnper == 0 ) pimplR->phvdper = ftmp1;
         else if( pimplR->phvnper == 1 ) pimplR->phvdper = ftmp1 - pimplR->phvdper;
         (pimplR->phvnper)++;
      }
      fin.close();
      pimplR->fphvnmem = fphvname;
   }

   // run rad_pattern
   if( type == 'R' ) {
      rad_pattern_r_( pimplR->feig_buff, &(pimplR->feig_len), &(pimplR->phvnper), &(pimplR->phvdper),
                      &(strike), &(dip), &(rake), &(depth), &(perlst.at(0)), &nper, azi, grT, phT, amp );
   } else if( type == 'L' ) {
      rad_pattern_l_( pimplR->feig_buff, &(pimplR->feig_len), &(pimplR->phvnper), &(pimplR->phvdper),
                      &(strike), &(dip), &(rake), &(depth), &(perlst.at(0)), &nper, azi, grT, phT, amp );
   }

   per_azi_pred.clear(); per_azi_pred.resize(nper);
   for(int i=0; i<nper; i++) {
      per_azi_pred.at(i).resize(181);
      for(int iazi=0; iazi<181; iazi++) {
         if( azi[iazi] != iazi*2 ) return false;
         per_azi_pred[i][iazi] = AziData( azi[iazi], grT[i][iazi], phT[i][iazi], amp[i][iazi] );
      }
   }

   return true;
}
