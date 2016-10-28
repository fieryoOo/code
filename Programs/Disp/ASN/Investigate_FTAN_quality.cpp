#include "DisAzi.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>

static float ghlgrv, ghlphv, alphagrv, alphaphv;

struct PointC {
	float x, y;
	PointC( const float x = NaN, const float y = NaN ) : x(x), y(y) {}
	PointC( const std::string &line, int ix=1, int iy=2 ) {
		if( ix >= iy )
			throw std::runtime_error("Error(PointC::PointC): ix must be <iy!");
		std::stringstream ss(line);
		for(int i=1; i<=iy; ++i) {
			if( ! (ss >> y) )
				throw std::runtime_error("Error(PointC::PointC): format error in line "+line);
			if( i == ix ) x = y;
		}
		//std::cerr<<x<<" "<<y<<std::endl;
	}
	friend bool operator< ( const PointC &p1, const PointC &p2 ) {
      return (p1.x < p2.x);
   }
	friend std::ostream &operator<<( std::ostream &o, const PointC &p) {
		o<<p.x<<" "<<p.y;
		return o;
	}
	static constexpr float NaN = -12345.;
};

struct STATION {
	std::string name;
	float lon, lat;
	STATION( const std::string &line ) {
		std::stringstream ss(line);
		if( ! ( ss >> name >> lon >> lat) )
         throw std::runtime_error("Error(STATION::STATION): format error in line "+line);
		if( lon < 0 ) lon += 360.;
		//std::cerr<<name<<" "<<lon<<" "<<lat<<std::endl;
	}
	friend std::ostream &operator<<( std::ostream &o, const STATION &sta) {
		o<<sta.name<<" "<<sta.lon<<" "<<sta.lat;
		return o;
	}
};

/*
struct LAG_DATA {
   float x;
   float disgrv, disphv;
   float amp, snr;
   float grv, phv;
	LAG_DATA( const std::string &line ) {
		throw std::runtime_error("Error(LAG_DATA::LAG:DATA): two line info needed!"); 
	}
	LAG_DATA( const std::string &displine, const std::string &snrline ) {
		float ftmp;
		if( sscanf(displine.c_str(), "%f %f %f %f %f", &ftmp, &ftmp, &x, &grv, &phv) != 5 ) 
			throw std::runtime_error("Error(LAG_DATA::LAG_DATA): format error in disp line "+displine);
      if( sscanf(snrline.c_str(), "%f %f %f", &ftmp, &amp, &snr) != 3 )
			throw std::runtime_error("Error(LAG_DATA::LAG_DATA): format error in snr line "+snrline);
		if( ftmp != x )
			throw std::runtime_error("Error(LAG_DATA::LAG_DATA): period mismatch between "+displine+" and "+snrline);
		x = log(x);
	}
	friend bool operator< ( const LAG_DATA &ld1, const LAG_DATA &ld2 ) {
      return (ld1.x < ld2.x);
   }
	friend std::ostream &operator<<( std::ostream &o, const LAG_DATA &d) {
		o<<d.x<<" "<<d.grv<<" "<<d.phv<<" "<<d.amp<<" "<<d.snr;
		return o;
	}
};
*/

template <class T>
class Curve {
public:
      // con/destructors
      // load from file
      Curve( const std::string& fname = "", int ix = 1, int iy = 2 ) {
         if( !fname.empty() ) Load(fname, ix, iy);
      }
      // initialized to a single value
      Curve( const size_t size, const T Tval = T{} ) {
         //const PointC& pc = T{};
         dataV.assign(size, Tval);
      }

		// load curve from file
		void Load( const std::string& fname, int ix, int iy ) {
			// check infile
			std::ifstream fin(fname);
			if( ! fin )
				throw std::runtime_error("Error(Curve::Load): IO failed on file "+fname);
			// read
			for(std::string line; std::getline(fin, line); ) {
				try {
					T t(line, ix, iy);
					dataV.push_back(t);
				} catch( std::exception& e ) {
					break;
				}
			}
			//std::cout<<dataV.size()<<" points loaded from file "<<fname<<std::endl;
			// sort
			sort();
		}

		void push_back(const T &t) { dataV.push_back(t); }

		void sort() { std::sort(dataV.begin(), dataV.end()); }

		float xfront() const { return dataV.front().x; }
		float xback() const { return dataV.back().x; }

		void clear() { dataV.clear(); }
      size_t size() const { return dataV.size(); }
      void reserve( size_t size ) { dataV.reserve(size); }

      typename std::vector<T>::const_iterator begin() const { return dataV.begin(); }
      typename std::vector<T>::const_iterator end() const { return dataV.end(); }
      typename std::vector<T>::iterator begin() { return dataV.begin(); }
      typename std::vector<T>::iterator end() { return dataV.end(); }

		float Val( const float x, bool allowOOB=false ) const {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), T(x) );
			// greater than end
			const float ferr = 1.0e-3;
			if( itupp >= dataV.end() ) {
				if( allowOOB || fabs((*(itupp-1)).x - x)<ferr ) return (*(itupp-1)).y;
				else return T::NaN;
			}
			// smaller than beg
			if( itupp<dataV.begin()+1 ) {
				if( allowOOB || fabs((*(itupp)).x - x)<ferr ) return (*(itupp)).y;
				return T::NaN;
			}
			const auto itlow = itupp - 1;
			if( itlow->x > x+ferr ) return T::NaN;
			// between two points: interpolate
			if( (*itlow).y!=T::NaN && (*itupp).y!=T::NaN ) {
				float vdiff = (*itupp).y - (*itlow).y, Tdiff = (*itupp).x - (*itlow).x;
				float deriv = vdiff / Tdiff;
				return (*itlow).y + ( x-(*itlow).x ) * deriv;
			}
			//if( (*itlow).y!=T::NaN ) return (*itlow).y;
			//if( (*itupp).y!=T::NaN ) return (*itupp).y;
			return T::NaN;
		}

		void cut( const float xl, const float xu, const int nmargin = 1 ) {
			sort();
			auto itlo = std::lower_bound( dataV.begin(), dataV.end(), T(xl) );
			dataV.erase( dataV.begin(), std::prev(itlo,nmargin) );
			auto itup = std::upper_bound( dataV.begin(), dataV.end(), T(xu) );
			dataV.erase( std::next(itup,nmargin), dataV.end() );
			//iter = dataV.begin();
		}

		Curve<T> DistS(const Curve<T> &c2) const {
			Curve<T> cDiff = *this;
			//std::cerr<<"debug: c1_range="<<xfront()<<"-"<<xback()<<"   c2_range="<<c2.xfront()<<"-"<<c2.xback()<<std::endl;
			for( auto &p : cDiff ) {
				float y2= c2.Val(p.x);
				if( y2 == T::NaN )
					throw std::runtime_error("Error(Curve::DistS): c2 does not cover c1!");
				float dismin = log(p.y/y2); dismin = dismin*dismin;
				for( auto &p2 : c2 ) {
					float dx = log(p2.x/p.x);
					float dy = log(p2.y/p.y);
					float dis = dx*dx+dy*dy;
					if( dis < dismin ) dismin = dis;
				}
				p.y = dismin;
			}
			return cDiff;
		}

		friend std::ostream &operator<<( std::ostream &o, const Curve<T> &c) {
			for( auto &t : c ) o << t <<"\n";
			return o;
		}

		// math operators
      Curve<T> operator-() const { return YTransform( std::negate<float>() ); }
      Curve<T> sqrt() const { return YTransform( [](const float val){ return std::sqrt(val); } ); }
      Curve<T> exp() const { return YTransform( [](const float val){ return std::exp(val); } ); }
      Curve<T> FillHoles( const float ymin ) const { return YTransform( [&](const float val){ return std::max(val, ymin); } ); }
      Curve<T>& operator+=( const Curve<T>& c2 ) { *this = (*this) + c2; return *this; }
      Curve<T>& operator-=( const Curve<T>& c2 ) { *this = (*this) - c2; return *this; }
      Curve<T>& operator*=( const Curve<T>& c2 ) { *this = (*this) * c2; return *this; }
      Curve<T>& operator/=( const Curve<T>& c2 ) { *this = (*this) / c2; return *this; }
      Curve<T>& operator*=( const float mul ) {
         for( auto& p : dataV ) p.y *= mul;
         return *this;
      }
      Curve<T>& operator/=( const float div ) {
         return ( (*this) *= (1./div) );
      }

      friend Curve<T> sqrt( const Curve<T>& c ) { return c.sqrt(); }
      friend Curve<T> exp( const Curve<T>& c ) { return c.exp(); }
      friend Curve<T> operator+(const Curve<T>& c1, const Curve<T>& c2) {
         return Arithmetic( c1, c2, std::plus<float>() );
      }
      friend Curve<T> operator-(const Curve<T>& c1, const Curve<T>& c2) {
         return Arithmetic( c1, c2, std::minus<float>() );
      }
      friend Curve<T> operator*(const Curve<T>& c1, const Curve<T>& c2) {
         return Arithmetic( c1, c2, std::multiplies<float>() );
      }
      friend Curve<T> operator*(const Curve<T>& c1, const float mul) {
         return c1.YTransform( [&](const float val){ return val*mul; } );
      }
      friend Curve<T> operator/(const Curve<T>& c1, const Curve<T>& c2) {
         return Arithmetic( c1, c2, std::divides<float>() );
      }
      friend Curve<T> operator/(const Curve<T>& c1, const float div) {
         return c1 * (1./div);
      }

		template<class Functor>
      Curve<T> YTransform( Functor func ) const {
		   Curve<T> c2 = *this;
		   //for_each( c2.dataV.begin(), c2.dataV.end(), func );
		   for( auto& p : c2.dataV ) 
				if(p.y!=T::NaN) p.y = func(p.y);
		   return c2;
		}

		template<class Functor>
      friend Curve<T> Arithmetic(const Curve<T>& c1, const Curve<T>& c2, Functor func) {
         Curve<T> c3;
         //std::cerr<<"In Arithmetic: c1size="<<c1.size()<<" c2size="<<c2.size()<<std::endl;
			/*
         if( c1.dataV.empty() ) {
            c3.reserve( c2.size() );
            for( auto p2 : c2.dataV ) {
               p2.y = func( 0., p2.y );
               c3.push_back(p2);
            }
         } else {
			*/
			c3.reserve( c1.size() );
			for( auto p1 : c1 ) {
				float val2 = c2.Val(p1.x);
				if( p1.y==T::NaN || val2==T::NaN ) {
					//throw ErrorCv::OutOfRange( FuncName, "x value = "+std::to_string(p1.x) );
					p1.y = T::NaN;
				} else {
					p1.y = func( p1.y, val2 );
				}
				c3.push_back(p1);
			}
			return c3;
		}

private:
		std::vector<T> dataV;
};


void AmpSNRSummation( const Curve<PointC> &Cgrv, const Curve<PointC> &Cphv, Curve<PointC> &Csnr, Curve<PointC> &Camp, const Curve<PointC> &CpredG, const Curve<PointC> &CpredP, float perl, float perh, float &snr, float &amp, const std::string &suf ) {

	//std::ofstream fpG("debug_predG.txt"); fpG<<CpredG;
	// compute shift/distance square in the log(per)-log(vel) space
	auto CDistSG = Cgrv.DistS(CpredG);	
	auto CDistSP = Cphv.DistS(CpredP);
	// compute weight from CDsitSG and CDistSP
	//auto CWeight = (CDistSG*alphagrv).exp() * 0.7 + (CDistSP*alphaphv).exp() * 0.3;
	auto CWeight = (CDistSG*alphagrv).exp();
	//std::ofstream fdebug1("GroupWeight_"+suf+".txt"); fdebug1<<(CDistSG*alphagrv).exp();

	// weight down amp and snr by CDist
	Camp = Camp.FillHoles(0.) * CWeight; 
	Csnr = Csnr.FillHoles(0.) * CWeight;
	//std::ofstream fdebug2("SNR_Downweighted_"+suf+".txt"); fdebug2<<Csnr;

	// define computing period range
	//if( ghlgrv<1. ) perl = std::max(perl, CpredG.xfront());
	//if( ghlgrv<1. ) perh = std::min(perh, CpredG.xback());

	// integrate amp and snr in [perl, perh]
	snr = amp = 0.;
	float xstep = 0.05, xb = std::log(perl)+xstep*0.5, x; 
	for(x=xb; x<=std::log(perh); x+=xstep) {
		float per = std::exp(x), snrcur = Csnr.Val(per), ampcur = Camp.Val(per);
		//std::cerr<<per<<" "<<x<<" "<<snrcur<<" "<<ampcur<<std::endl;
		if(snrcur>0.) snr += snrcur * xstep;
		if(ampcur>0.) amp += ampcur * xstep;
	}
	snr /= (x-xb); amp /= (x-xb);
}

float TimeShift( const Curve<PointC> &Cv1, const Curve<PointC> &Cv2, Curve<PointC> &Csnr1, Curve<PointC> &Csnr2, 
		const float dist, const float perl, const float perh, const float snrmin ) {
	//auto TDiff = (dist/Cv1) - (dist/Cv2);
	auto TDiff = (Cv2-Cv1) / (Cv2*Cv1) * dist;
	// integrate TDiff in [perl, perh]
	int neff = 0;
	float xstep = 0.05, Tshift = 0.; 
	for(float x=std::log(perl); x<=std::log(perh); x+=xstep) {
		float per = std::exp(x), snr1 = Csnr1.Val(per), snr2 = Csnr2.Val(per);
		//std::cerr<<per<<" "<<x<<" "<<TDiff.Val(per)<<" "<<snr1<<" "<<snr2<<std::endl;
		if( snr1>snrmin && snr2>snrmin ) {
			Tshift += TDiff.Val(per); ++neff;
		}
	}
	return neff==0 ? PointC::NaN : (Tshift / neff);
}


int main(int argc, char *argv[]) {
	if(argc!=9) {
		std::cerr<<"usage: "<<argv[0]<<" [station.lst] [stanm1-stanm2-fDisppos-fAmppos-fDispneg-fAmpneg-daynum-fpredDispgrv-fpredDispphv list] [period_lowerend] [period_higherend] [ghlgrv (ignore predictions if >=1)] [ghlphv (ignore predictions if >=1)] [SNRmin(only for Tshift)] [outname]"<<std::endl;
		return -1;
	}

	// get parameters
	std::string stalst(argv[1]), pathlst(argv[2]), outname(argv[8]);
	ghlgrv = atof(argv[5]); ghlphv = atof(argv[6]);
	alphagrv = -0.5/(ghlgrv*ghlgrv); alphaphv = -0.5/(ghlphv*ghlphv);
	float perl = atof(argv[3]), perh = atof(argv[4]);
	float snrmin = atof(argv[7]); //, wsnr_pow = atoi(argv[8]);
	if( ghlgrv >= 1. ) std::cout<<"*** Warning(main): Group Predictions will be ignored! ***"<<std::endl;
	if( ghlphv >= 1. ) std::cout<<"*** Warning(main): Phase Predictions will be ignored! ***"<<std::endl;
	if( snrmin < 0. ) {
		snrmin = 0.;
		std::cout<<"*** Warning(main): snrmin cannot be negative, corrected back to 0 ***"<<std::endl;
	}

	// load in the station list
	std::unordered_map<std::string, STATION> staM;
	std::ifstream fsta(stalst);
	if( ! fsta ) throw std::runtime_error("Error(main): IO failed on file "+stalst);
	for(std::string line; std::getline(fsta, line); ) {
		STATION sta(line);
		staM.emplace(sta.name, sta);
	}
	std::cout<<"### "<<staM.size()<<" stations loaded. ###"<<std::endl;

	// output file
	std::ofstream fout(outname);
	fout<<"sta1(1) lon1(2) lat1(3)  sta2(4) lon2(5) lat2(6)  dist(7) azi1(8) azi2(9) dnum(10) : snrsig_pos(12) ampsig_pos/dnum(13) snrsig_neg(14) ampsig_neg/dnum(15) GrvTimeShift(16) PhvTimeShift(17)\n";

	// read the path list by line
	std::ifstream fpath(pathlst);
	if( ! fpath ) throw std::runtime_error("Error(main): IO failed on file "+pathlst);
	for(std::string line; std::getline(fpath,line); ) {
		// phrase the input line
		std::string stanm1, stanm2, disp_pf, snr_pf, disp_nf, snr_nf, pdisp_g, pdisp_p;
		int daynum; std::stringstream ss(line);
		if( ! (ss>>stanm1>>stanm2>>disp_pf>>snr_pf>>disp_nf>>snr_nf>>daynum>>pdisp_g>>pdisp_p) ) {
			std::cerr<<"Warning(main): format error in line "<<line<<" from file "<<pathlst<<std::endl;
			continue;
		}
		//std::cerr<<stanm1<<" "<<stanm2<<" "<<disp_pf<<" "<<snr_pf<<" "<<disp_nf<<" "<<snr_nf<<" "<<daynum<<" "<<pdisp_g<<" "<<pdisp_p<<std::endl;
		// start working on the current path
		try {
			// get station info from staM
			STATION &sta1 = staM.at(stanm1), &sta2 = staM.at(stanm2);
			std::cout<<"### Extracting information for path "<<stanm1<<" - "<<stanm2<<std::endl;
			Path<float> path(sta1.lon, sta1.lat, sta2.lon, sta2.lat);
			// load curves
			Curve<PointC> predG(pdisp_g), predP(pdisp_p);
			Curve<PointC> grvp(disp_pf, 3, 4), phvp(disp_pf, 3, 5);
			Curve<PointC> ampp(snr_pf, 1, 2), snrp(snr_pf, 1, 3);
			Curve<PointC> grvn(disp_nf, 3, 4), phvn(disp_nf, 3, 5);
			Curve<PointC> ampn(snr_nf, 1, 2), snrn(snr_nf, 1, 3);
			// trim
			float pl = perl, pu = perh;
			grvp.cut(pl,pu,2); phvp.cut(pl,pu,2);	ampp.cut(pl,pu,2); snrp.cut(pl,pu,2);
			grvn.cut(pl,pu,2); phvn.cut(pl,pu,2);	ampn.cut(pl,pu,2); snrn.cut(pl,pu,2);
			// correct 2pi for phase
			float dist = path.Dist();
			for( auto &p : phvp )
				p.y = dist / (dist/p.y - floor((dist/p.y-dist/predP.Val(p.x))/p.x+0.5) * p.x);
			for( auto &p : phvn )
				p.y = dist / (dist/p.y - floor((dist/p.y-dist/predP.Val(p.x))/p.x+0.5) * p.x);
			// compute the effective(average-down-weighted) amp and snr
			float snrsig_pos, snrsig_neg, ampsig_pos, ampsig_neg;
			AmpSNRSummation( grvp, phvp, snrp, ampp, predG, predP, perl, perh, snrsig_pos, ampsig_pos, "pos" );
			AmpSNRSummation( grvn, phvn, snrn, ampn, predG, predP, perl, perh, snrsig_neg, ampsig_neg, "neg" );
			// compute time shifts on segments with an effective snr>snrmin
			// note that the effective snrs are already downweighted by diff(grv,grvpred), so the resulting Tshifts are averaged only on segs with good data
			float TshiftG = TimeShift( grvp, grvn, snrp, snrn, dist, perl, perh, snrmin );
			float TshiftP = TimeShift( phvp, phvn, snrp, snrn, dist, perl, perh, snrmin );
			fout<<sta1<<"   "<<sta2<<"   "<<dist<<" "<<path.Azi1()<<" "<<path.Azi2()<<" "<<daynum<<" : "<<snrsig_pos<<" "<<ampsig_pos/daynum<<"   "<<snrsig_neg<<" "<<ampsig_neg/daynum<<"   "<<TshiftG<<" "<<TshiftP<<std::endl;
		} catch(std::exception &e){
			std::cerr<<"Exception caught: "<<e.what()<<std::endl;
		}
	}

	return 0;
}

