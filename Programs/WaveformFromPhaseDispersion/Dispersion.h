#ifndef DISPERSION_H
#define DISPERSION_H

#include "Curve.h"

class Dispersion;
class Ddata : public PointC {
public:
	Ddata( float xin=NaN, float yin=NaN, float sdensin=1. )
		: PointC(xin, yin, sdensin) {
		if(x != NaN ) {
			om = twopi / x;
			if( y != NaN ) k_ang = om / y;
		}
	}

	Ddata(const std::string& input, int ix=1, int iy=2, int iz=3) 
		: PointC(input, ix, iy, iz) {
		om = twopi / x;
		k_ang = om / y;
	}

	friend Curve<Ddata>;
	friend Dispersion;

private:
	float om = NaN, k_ang = NaN;	//angular wavenumber
};

// to save space, use PointC.z to store om
class KDeriv : public Curve<PointC> {
public:
		float Deriv_2om( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), PointC(per) );
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return PointC::NaN;
			const auto itlow = itupp - 1;
			return ( (*itupp).y - (*itlow).y ) / ( (*itupp).z - (*itlow).z );
		}

		//float Reciprocal( std::vector<PointC>& grvV ) {
		void Reciprocal( KDeriv& grvs ) {
			//grvV.resize( dataV.size() );
			grvs.clear(); grvs.reserve( dataV.size() );
			for(int i=0; i<dataV.size(); i++) {
				grvs.push_back(PointC(dataV[i].x, 1./dataV[i].y, dataV[i].z));
				//grvV[i] = PointC(dataV[i].x, 1./dataV[i].vel, dataV[i].z);
			}
		}
};

class Dispersion : public Curve<Ddata> {
public:
		Dispersion()
			: Curve() {}

		Dispersion( const std::string& fname )
			: Curve(fname) {}

		void push_back( float perin, float velin, float sdensin=1. ) {
			dataV.push_back( Ddata(perin, velin, sdensin) );
			//std::sort(dataV.begin(), dataV.end());
		}

		//void Sort() { std::sort(dataV.begin(), dataV.end()); }

		float Om( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), Ddata(per) );
			if( itupp == dataV.end() ) {
				if( fabs((*(itupp-1)).x - per) < 0.1 ) return (*(itupp-1)).om;
				else return Ddata::NaN;
			}
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return Ddata::NaN;
			const auto itlow = itupp - 1;
			float odiff = (*itupp).om - (*itlow).om, Tdiff = (*itupp).x - (*itlow).x;
			float deriv = odiff / Tdiff;
			return (*itlow).om + ( per-(*itlow).x ) * deriv;
		}

		float Sdens( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), Ddata(per) );
			if( itupp == dataV.end() ) {
				if( fabs((*(itupp-1)).x - per) < 0.1 ) return (*(itupp-1)).z;
				else return Ddata::NaN;
			}
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return Ddata::NaN;
			const auto itlow = itupp - 1;
			float sdiff = (*itupp).z - (*itlow).z, Tdiff = (*itupp).x - (*itlow).x;
			float deriv = sdiff / Tdiff;
			return (*itlow).z + ( per-(*itlow).x ) * deriv;
		}

		float Deriv_k2om( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), Ddata(per) );
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return Ddata::NaN;
			const auto itlow = itupp - 1;
			return ( (*itupp).k_ang - (*itlow).k_ang ) / ( (*itupp).om - (*itlow).om );
		}

		bool Deriv_k2om( KDeriv& curveout ) {
			curveout.clear();
			curveout.reserve( dataV.size() );
			if( dataV.size() <= 1 ) return false;
			for( auto it=dataV.begin()+1; it<dataV.end(); it++ ) {
				float per = 0.5 * ( (*it).x + (*(it-1)).x );
				float om = Om( per );
				float deriv = Deriv_k2om( per );
				curveout.push_back(PointC(per, deriv, om));
			}
		}

		void ComputeSpectrum( const std::string& sacname, float tmin, float tmax ) {
			SacRec sac(sacname);
			sac.Load(); sac.cut(tmin, tmax);
			SacRec sac_am;
			sac.ToAm(sac_am);
			sac_am.Smooth(0.0002, sac);
			float *sigsac = sac.sig.get();
			for( auto& data : dataV ) {
				float freq = 1./data.x;
				int ifreq = (int)floor(freq/sac.shd.delta+0.5);
				data.z = sigsac[ifreq];
				//std::cout<<freq<<" "<<data.z<<std::endl;
			}
		}

};

#endif
