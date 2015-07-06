#ifndef MODELINFO_H
#define MODELINFO_H

#include "MyOMP.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <random>
#include <chrono>
#include <thread>

/* -------------------- the RNG class-------------------- */
class Rand {
   std::default_random_engine generator1;
   std::uniform_real_distribution<float> d_uniform;
   std::normal_distribution<float> d_normal;
public:
   Rand() /* add a true random number from std::random_device to the time seed to ensure thread safety */
      : generator1( std::chrono::system_clock::now().time_since_epoch().count() + std::random_device{}() )
      , d_uniform(0., 1.)
      , d_normal(0., 1.) {}

   //~Rand() {} // should not be defined!!!

   inline float Uniform() { return d_uniform(generator1); }
   inline float Normal() { return d_normal(generator1); }

};

template< class T >
struct FocalInfo {
	static constexpr float NaN = -12345.;
   T stk, dip, rak, dep;

   //FocalInfo( T stkin = 180, T dipin = 45, T rakin = 0, T depin = 10 )
   FocalInfo( T stkin = NaN, T dipin = NaN, T rakin = NaN, T depin = NaN )
      : stk(stkin), dip(dipin), rak(rakin), dep(depin) {}

	/* check validation */
	virtual bool isValid() const {
		return (stk>=0.&&stk<360.) && (dip>=0.&&dip<=90.) &&
				 (rak>=-180.&&rak<180.) && (dep>=0.) ;
	}

	/* correct into the right range */
	inline void CorrectRange() {
		Correct2PI();
		if( stk==NaN || rak==NaN || dip==NaN || dip>=180. || dip<=-90. )
			throw std::runtime_error("Error(FocalInfo::CorrectRange): Invalid finfo!");
		if( dip >= 90. ) dip = 180. - dip;
		if( dip < 0. ) dip = -dip;
	}

	/* Transfer to the auxiliary nodal plane and slip of the current */
	void Auxiliary() {
		float s2, d2, r2;
		// convert deg to rad
		float deg2rad = M_PI/180.;
		stk *= deg2rad; dip *= deg2rad; rak *= deg2rad;
		// dip
		d2 = acos( sin(rak)*sin(dip) );
		// rak
		float sin_r2 = cos(dip) / sin(d2);
		float sin_dphi = cos(rak) / sin(d2);
		float cos_r2 = -sin(dip) * sin_dphi;	// 0 ~ pi
		r2 = asin(sin_r2);	// -pi/2 ~ pi/2
		if( cos_r2 < 0. ) r2 = M_PI - r2;
		// stk
		float cos_dphi = -1. / (tan(dip)*tan(d2));
		float dphi = asin(sin_dphi);
		if( cos_dphi < 0. ) dphi = M_PI - dphi;
		s2 = stk - dphi;
		// check and correct quadrant
		if( d2 > 0.5*M_PI ) {
			s2 += M_PI;
			d2 = M_PI - d2;
			r2 = 2.*M_PI - r2;
		}
		// convert back;
		float rad2deg = 1./deg2rad;
		stk = s2*rad2deg; dip = d2*rad2deg; rak = r2*rad2deg;
		Correct2PI();
	}

	void MomentTensor( const float M0 ) const {
		float deg2rad = M_PI/180.;
		float stk = this->stk, dip = this->dip, rak = this->rak;
		stk *= deg2rad; dip *= deg2rad; rak *= deg2rad;
		float sins = sin(stk), coss = cos(stk), sin2s = sin(2.*stk), cos2s = cos(2.*stk);
		float sind = sin(dip), cosd = cos(dip), sin2d = sin(2.*dip), cos2d = cos(2.*dip);
		float sinr = sin(rak), cosr = cos(rak);
		const float Mxx = -M0 * (sind*cosr*sin2s + sin2d*sinr*sins*sins);
		const float Mxy =  M0 * (sind*cosr*cos2s + sin2d*sinr*sins*coss);
		const float Mxz = -M0 * (cosd*cosr*coss + cos2d*sinr*sins);
		const float Myy =  M0 * (sind*cosr*sin2s - sin2d*sinr*coss*coss);
		const float Myz = -M0 * (cosd*cosr*sins - cos2d*sinr*coss);
		const float Mzz =  M0 * (sin2d*sinr);
		std::cout<<Mxx<<" "<<Mxy<<" "<<Mxz<<" "<<Myy<<" "<<Myz<<" "<<Mzz<<std::endl;
	}

   friend std::ostream& operator<< ( std::ostream& o, const FocalInfo& f ) {
      //o<<std::fixed<<std::setprecision(2)<<f.stk<<" "<<std::setw(6)<<f.dip<<" "<<std::setw(6)<<f.rak<<"  "<<std::setw(6)<<f.dep;
      o<<std::fixed<<std::setprecision(3)
		 <<std::setw(7)<<f.stk<<" "<<std::setw(7)<<f.dip<<" "<<std::setw(8)<<f.rak<<" "<<std::setw(7)<<f.dep; 
      return o; 
   }

   friend bool operator== ( const FocalInfo<T>& fi1, const FocalInfo<T>& fi2 ) {
      T dis_st = fabs(fi1.stk - fi2.stk);
      T dis_di = fabs(fi1.dip - fi2.dip);
      T dis_ra = fabs(fi1.rak - fi2.rak);
      T dis_de = fabs(fi1.dep - fi2.dep);
      return (dis_st<0.1 && dis_di<0.1 && dis_ra<0.1 && dis_de<0.1);
   }

private:
	/* shift by 2PIs into the correct range */
	void Correct2PI() {
		if( stk==NaN || rak==NaN ) return;
		while( stk >= 360. ) stk -= 360.;
		while( stk < 0. ) stk += 360.;
		while( rak >= 180. ) rak -= 360.;
		while( rak < -180. ) rak += 360.;
	}

};
typedef float ftype;

/* earthquake epicenter information */
struct EpicInfo {
	static constexpr float NaN = FocalInfo<ftype>::NaN;
	float lon, lat, t0;

	EpicInfo( float lonin = NaN, float latin = NaN, float t0in = NaN )
		: lon(lonin), lat(latin), t0(t0in) {}

	virtual bool isValid() const {
		return ( (lon<360.&&lon>=0.) && (lat>=-90.&&lat<=90.) && t0!=NaN );
	}

	friend std::ostream& operator<< ( std::ostream& o, const EpicInfo& e ) {
		o<<std::fixed<<std::setprecision(4)<<e.lon<<" "<<e.lat<<"  "<<std::setw(7)<<e.t0; 
		return o; 
	}

	friend bool operator== ( const EpicInfo& ei1, const EpicInfo& ei2 ) {
		float dis_lon = fabs(ei1.lon - ei2.lon) * 100.;
		float dis_lat = fabs(ei1.lat - ei2.lat) * 100.;
		float dis_t = fabs(ei1.t0 - ei2.t0);
		return (dis_lon<0.01 && dis_lat<0.01 && dis_t<0.01);
	}
};


/* model space */
class ModelInfo : public FocalInfo<ftype>, public EpicInfo {
	public:
		using FocalInfo::NaN;

		ModelInfo() {}

		ModelInfo( const float lonin, const float latin, const float timin,
				const float stkin, const float dipin, const float rakin, const float depin )
			: EpicInfo(lonin, latin, timin)
			, FocalInfo(stkin, dipin, rakin, depin) {}

		ModelInfo( const EpicInfo& einfo, const FocalInfo& finfo )
			: EpicInfo(einfo), FocalInfo(finfo) {}

		virtual bool isValid() const {
			return ( FocalInfo<ftype>::isValid() && EpicInfo::isValid() );
		}

		friend std::ostream& operator<< ( std::ostream& o, const ModelInfo& m ) {
			o<< static_cast< const FocalInfo<ftype>& >(m) << "   " << static_cast< const EpicInfo& >(m);
			return o;
		}

		friend bool operator== ( const ModelInfo& ms1, const ModelInfo& ms2 ) {
			return ( ( static_cast< const FocalInfo<ftype>& >(ms1) == static_cast< const FocalInfo<ftype>& >(ms2) )
					&& ( static_cast< const EpicInfo& >(ms1) == static_cast< const EpicInfo& >(ms2) ) );
		}
};

#endif
