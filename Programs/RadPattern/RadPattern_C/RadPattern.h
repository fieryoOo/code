#ifndef RADPATTERN_H
#define RADPATTERN_H

#include "MyOMP.h"
#include <cmath>
#include <memory>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>


/* ---------- exceptions ---------- */
//#define FuncName __PRETTY_FUNCTION__
#define FuncName __FUNCTION__
namespace ErrorRP {
   class BadFile : public std::runtime_error {
   public:
      BadFile(const std::string funcname, const std::string info = "")
	 : runtime_error("Error("+funcname+"): Cannot access file ("+info+").") {}
   };

   class BadParam : public std::runtime_error {
   public:
      BadParam(const std::string funcname, const std::string info = "")
        : runtime_error("Error("+funcname+"): Bad parameters ("+info+").") {}
   };

   class BadAzi : public std::runtime_error {
   public:
      BadAzi(const std::string funcname, const std::string info = "")
        : runtime_error("Error("+funcname+"): Unexpected azimuths ("+info+").") {}
   };

   class BadBuff : public std::runtime_error {
   public:
      BadBuff(const std::string funcname, const std::string info = "")
        : runtime_error("Error("+funcname+"): Internal buffer modified ("+info+"). The RadPattern class is not yet thread-safe!") {}
   };

};


/*
struct AziData {
   bool valid;
   float azi;
   float misG, varG;
   float misP, varP;
   float A, varA;
   AziData()
      : azi(-12345.), misG(-12345.), misP(-12345.), A(-12345.) 
      , valid(false), varG(-12345.), varP(-12345.), varA(-12345.) {}

   AziData( float aziin, float misGin, float misPin, float ampin )
      : azi(aziin), misG(misGin), misP(misPin), A(ampin)
      , valid(true), varG(-12345.), varP(-12345.), varA(-12345.) {}

	bool isGood() const { return ( valid && misG!=NaN && varG!=NaN && misP!=NaN && varP!=NaN && A!=NaN ); }

	friend bool operator<( const AziData& ad1, const AziData& ad2 ) {
      return (ad1.azi < ad2.azi);
   }

   friend std::ostream& operator<< ( std::ostream& o, const AziData& ad ) {
      o<<"( "<<ad.azi<<"  "<<ad.misG<<" "<<ad.varG<<"  "<<ad.misP<<" "<<ad.varP
       <<"  "<<ad.A<<" "<<ad.varA<<" )";
      return o;
   }

protected:
	static constexpr float NaN = -12345.;
};
*/

typedef float ftype;

class RadPattern {
public:
	char type;
	float stk, dip, rak, dep, M0;

   RadPattern();
   RadPattern( const RadPattern& );
   RadPattern( RadPattern&& );
   RadPattern& operator= ( const RadPattern& );
   RadPattern& operator= ( RadPattern&& );
   ~RadPattern();

   /* Predict Rayleigh/Love wave radiation patterns */
   bool Predict( char type, const std::string& feigname, const std::string& fphvname,
					  const ftype strike, const ftype dip, const ftype rake,
					  const ftype depth, const ftype M0, const std::vector<float>& perlst );
				//std::vector< std::vector<AziData> >& per_azi_pred );

	/* get I0 at per */
	float I0( const float per ) const;

	/* prediction at one single azimuth. return false if the given azimuth is invalidated due to small amplitude */
	// M0 = scalar seismic momentum
	// dis = distance; alpha = attenuation coeff
	// J = mode energy integration (from eigen);
	// U = local group velocity
	bool GetPred( const float per, const float azi,	float& grt, float& pht, float& amp,
					  const float dis = NaN, const float alpha = NaN, const float J = NaN, const float U = NaN ) const;

	void OutputPreds( const std::string& fname, const float Afactor = 1. );

public:
	static constexpr float NaN = -12345.;
	static constexpr float oofourpi = 0.25e-13 / M_PI;	// unit convertion = 1.0e-13
	static constexpr int nazi = 181;
	static constexpr int dazi = 2;
	static constexpr int InvalidateHwidth = 3;	// half width in iazi of the focal pred invalidating window
	static constexpr float AmpValidPerc = 0.05;  // focal predictions with A < Aaverage*AmpValidPerc are invalidated

private:
   struct Rimpl;
   std::unique_ptr<Rimpl> pimplR;

	// sample azimuths
	std::vector<float> aziV;
	// I0 (mod energy integral) keyed by period
	std::map< float, float > I0M;
	// group, phase, amplitudes keyed by period
	std::map< float, std::vector<float> > grtM, phtM, ampM;

	void ShiftCopy( std::vector<float>& Vout, const float* arrayin, const int nazi ) const;
};

#endif
