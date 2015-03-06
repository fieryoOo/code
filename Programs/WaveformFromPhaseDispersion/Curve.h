#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

/* ----- single point data structures ----- */
template <class T> class Curve;
class KDeriv;
class Parabola;
class Point {
public:
	Point( float perin=NaN, float velin=NaN, float sdensin=1. ) 
		: per(perin), vel(velin), sdensity(sdensin) {}

	Point( const std::string& input ) {
		int nrd = sscanf(input.c_str(), "%f %f %f", &per, &vel, &sdensity);
		if( nrd < 2 )
			throw std::runtime_error("Error(Point::Point): format error");
	}

	friend bool operator< ( const Point& p1, const Point& p2 ) {
		return (p1.per < p2.per);
	}

	friend std::ostream& operator<< ( std::ostream& o, const Point& pt ) {
		o << pt.per << " " << pt.vel << " " << pt.sdensity;
		return o;
	}

	friend Curve<Point>;
	friend KDeriv;
	friend Parabola;
	friend Curve<Point> operator-(const Curve<Point>& c1, const Curve<Point>& c2);

	static constexpr float NaN = -12345.;
	static constexpr float twopi = M_PI * 2.;

//protected:
	float per = NaN, vel = NaN, sdensity = 1.;
};

class Dispersion;
class Ddata : public Point {
public:
	Ddata( float perin=NaN, float velin=NaN, float sdensin=1. )
		: Point(perin, velin, sdensin) {
		if(per != NaN ) {
			om = twopi / per;
			if( vel != NaN ) k_ang = om / vel;
		}
	}

	Ddata(const std::string& input) 
		: Point(input) {
		om = twopi / per;
		k_ang = om / vel;
	}

	friend Curve<Ddata>;
	friend Dispersion;

private:
	float om = NaN, k_ang = NaN;	//angular wavenumber
};

/* ----- data containers ----- */
template <class T>
class Curve {
public:
		Curve() {}
		Curve( const std::string& fname ) { Load(fname); }
		void Load( const std::string& fname ) {
			// check infile
			std::ifstream fin(fname);
			if( ! fin )
				throw std::runtime_error("Error(Curve::Load): cannot access dispersion infile " + fname);
			// read
			for(std::string line; std::getline(fin, line); ) {
				dataV.push_back( T(line) );
			}
			std::cout<<dataV.size()<<" points loaded from file "<<fname<<std::endl;
			// sort
			std::sort(dataV.begin(), dataV.end());
			// debug
			//for( const auto& d : dataV ) std::cerr<<d<<std::endl;
		}

		void clear() { dataV.clear();	}
		void push_back( float per, float vel, float om ){ dataV.push_back( Point(per, vel, om) ); }

		void rewind() { iter = dataV.begin(); }
		bool get( T& Tout ) {
			if( iter >= dataV.end() ) return false;
			Tout = *iter; return true;
		}
		void next() { iter++; }

		float xmin() { return dataV.front().per; }
		float xmax() { return dataV.back().per; }

		float Val( const float per ) const {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), T(per) );
			if( itupp == dataV.end() ) {
				if( fabs((*(itupp-1)).per - per) < 0.1 ) return (*(itupp-1)).vel;
				else return Ddata::NaN;
			}
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return T::NaN;
			const auto itlow = itupp - 1;
			float vdiff = (*itupp).vel - (*itlow).vel, Tdiff = (*itupp).per - (*itlow).per;
			float deriv = vdiff / Tdiff;
			return (*itlow).vel + ( per-(*itlow).per ) * deriv;
		}

		float Deriv( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), T(per) );
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return T::NaN;
			const auto itlow = itupp - 1;
			return ( (*itupp).vel - (*itlow).vel ) / ( (*itupp).per - (*itlow).per );
		}

		void FindValue(float val, std::vector<float>& perV) {
			perV.clear();
			if( dataV.size() <= 1 ) return;
			for( int i=0; i<dataV.size()-1; i++ ) {
				const auto& datacur = dataV[i];
				float velh = dataV[i+1].vel, vell = datacur.vel;
				if( val == vell ) {	// equal to or
					perV.push_back(datacur.per);
				} else if( (val-vell) * (val-velh) < 0. ) { // in between
					float perh = dataV[i+1].per, perl = datacur.per;
					float vdiff = velh - vell, pdiff = perh - perl;
					//std::cerr<<"looking for "<<val<<" in between "<<perl<<"-"<<perh<<"sec "<<vell<<"-"<<velh<<"s/km"<<std::endl;
					perV.push_back( perl + (val-vell) * pdiff / vdiff );
				}
			}
			const auto& datacur = dataV.back();
			if( val == datacur.vel ) perV.push_back(datacur.per);
		}

		void Write( const std::string& fname ) {
			std::ofstream fout( fname );
			if( ! fout ) {
				std::cerr<<"Error(Curve::Write): cannot write to file "<<fname<<std::endl;;
				exit(-2);
			}
			for( const auto& spP : dataV )
				fout<<spP<<"\n";
		}

		friend Curve<T> operator-(const Curve<T>& c1, const Curve<T>& c2) {
			Curve<T> c3;
			for( auto p1 : c1.dataV ) {
				float val2 = c2.Val(p1.per);
				if( val2 == Point::NaN ) continue;
				p1.vel -= val2;
				c3.dataV.push_back(p1);
			}
			return c3;
		}

protected:
		std::vector<T> dataV;
private:
		typename std::vector<T>::iterator iter = dataV.begin();
};

// to save space, use Point.sdensity to store om
class KDeriv : public Curve<Point> {
public:
		float Deriv_2om( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), Point(per) );
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return Point::NaN;
			const auto itlow = itupp - 1;
			return ( (*itupp).vel - (*itlow).vel ) / ( (*itupp).sdensity - (*itlow).sdensity );
		}

		//float Reciprocal( std::vector<Point>& grvV ) {
		float Reciprocal( KDeriv& grvs ) {
			//grvV.resize( dataV.size() );
			grvs.clear();
			for(int i=0; i<dataV.size(); i++) {
				grvs.push_back(dataV[i].per, 1./dataV[i].vel, dataV[i].sdensity);
				//grvV[i] = Point(dataV[i].per, 1./dataV[i].vel, dataV[i].sdensity);
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

		void Sort() { std::sort(dataV.begin(), dataV.end()); }

		float Om( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), Ddata(per) );
			if( itupp == dataV.end() ) {
				if( fabs((*(itupp-1)).per - per) < 0.1 ) return (*(itupp-1)).om;
				else return Ddata::NaN;
			}
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return Ddata::NaN;
			const auto itlow = itupp - 1;
			float odiff = (*itupp).om - (*itlow).om, Tdiff = (*itupp).per - (*itlow).per;
			float deriv = odiff / Tdiff;
			return (*itlow).om + ( per-(*itlow).per ) * deriv;
		}

		float Sdens( const float per ) {
			const auto itupp = std::upper_bound( dataV.begin(), dataV.end(), Ddata(per) );
			if( itupp == dataV.end() ) {
				if( fabs((*(itupp-1)).per - per) < 0.1 ) return (*(itupp-1)).sdensity;
				else return Ddata::NaN;
			}
			if( itupp<dataV.begin()+1 || itupp>=dataV.end() )
				return Ddata::NaN;
			const auto itlow = itupp - 1;
			float sdiff = (*itupp).sdensity - (*itlow).sdensity, Tdiff = (*itupp).per - (*itlow).per;
			float deriv = sdiff / Tdiff;
			return (*itlow).sdensity + ( per-(*itlow).per ) * deriv;
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
			if( dataV.size() <= 1 ) return false;
			for( auto it=dataV.begin()+1; it<dataV.end(); it++ ) {
				float per = 0.5 * ( (*it).per + (*(it-1)).per );
				float om = Om( per );
				float deriv = Deriv_k2om( per );
				curveout.push_back(per, deriv, om);
			}
		}

		void ComputeSpectrum( const std::string& sacname ) {
			SacRec sac(sacname);
			sac.Load();
			SacRec sac_am;
			sac.ToAm(sac_am);
			sac_am.Smooth(0.0002, sac);
			float *sigsac = sac.sig.get();
			for( auto& data : dataV ) {
				float freq = 1./data.per;
				int ifreq = (int)floor(freq/sac.shd.delta+0.5);
				data.sdensity = sigsac[ifreq];
				//std::cout<<freq<<" "<<data.sdensity<<std::endl;
			}
		}

};


class Parabola {
public:
	Parabola( const Point& P1in, const Point& P2in, const Point& P3in )
		: P1(P1in), P2(P2in), P3(P3in) {}

	void Solve(){
		float x1 = P1.per, y1 = P1.vel;
		float x2 = P2.per, y2 = P2.vel;
		float x3 = P3.per, y3 = P3.vel;
		float xs1 = x1*x1, xs2 = x2*x2, xs3 = x3*x3;
		float denom = (x1-x2) * (x1-x3) * (x2-x3);
		a = (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)) / denom;
		b = (xs3*(y1-y2) + xs2*(y3-y1) + xs1*(y2-y3)) / denom;
		c = (x2*x3*(x2-x3)*y1 + x3*x1*(x3-x1)*y2 + x1*x2*(x1-x2)*y3) / denom;
	}

	Point Vertex() {
		if(a==NaN || b==NaN || c==NaN)
			Solve();
		if( PV.per==NaN || PV.vel==NaN ) {
			PV.per = - b / (2.*a);
			PV.vel = c - a*PV.per*PV.vel;
		}
		return PV;
	}

	float A() { 
		if(a==NaN) Solve(); 
		return a;
	}
	float B() { 
		if(b==NaN) Solve(); 
		return b;
	}
	float C() { 
		if(c==NaN) Solve(); 
		return c;
	}

private:
	Point P1, P2, P3;
	Point PV;
	float NaN = Point::NaN;
	float a=NaN, b=NaN, c=NaN;
};
