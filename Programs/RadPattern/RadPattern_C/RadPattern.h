#ifndef RADPATTERN_H
#define RADPATTERN_H

#include <memory>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

struct FocalInfo {
   int strike, dip, rake;
   int depth;
   FocalInfo( int strikein = 180, int dipin = 45, int rakein = 0, int depthin = 10 )
      : strike(strikein), dip(dipin), rake(rakein), depth(depthin) {}
   friend std::ostream& operator<< ( std::ostream& o, const FocalInfo& f ) {
      o<<"( "<<std::setw(3)<<f.strike<<" "<<std::setw(2)<<f.dip<<" "<<std::setw(4)<<f.rake<<"  "<<std::setw(2)<<f.depth<<" )"; 
      return o; 
   }
};


struct AziData {
   float azi;
   float misG, varG;
   float misP, varP;
   float A, varA;
   AziData( float aziin = -12345., float misGin = -12345.,
            float misPin = -12345., float ampin = -12345. )
      : azi(aziin), misG(misGin), misP(misPin), A(ampin) {}

   friend std::ostream& operator<< ( std::ostream& o, const AziData& ad ) {
      o<<"( "<<ad.azi<<"  "<<ad.misG<<" "<<ad.varG<<"  "<<ad.misP<<" "<<ad.varP
       <<"  "<<ad.A<<" "<<ad.varA<<" )";
      return o;
   }

};


class RadPattern {
   struct Rimpl;
   std::unique_ptr<Rimpl> pimplR;

public:
   RadPattern();
   ~RadPattern();

   /* Predict Rayleigh wave radiation patterns */
   bool Predict( char type, const std::string& feigname, const std::string& fphvname, const FocalInfo& finfo,
                 std::vector<float>& perlst, std::vector< std::vector<AziData> >& per_azi_pred );

};

#endif
