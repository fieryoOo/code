#ifndef DISAZI_H
#define DISAZI_H

#include "Point.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
//using namespace std;

// in case Point.h not available
#ifndef POINT_H
#define POINT_H
template <class T>
class Point {
protected:
   T lon, lat;
public:
   Point(T lonin = -12345., T latin = -12345.)
      : lon(lonin), lat(latin) {}
	Point( const std::string& line ) {
		LoadLine(line);
	}
	bool LoadLine( const std::string& line ) {
		return ( sscanf(line.c_str(), "%f %f", &lon, &lat) == 2 );
	}

   inline const T& Lat() const { return lat; }
   inline       T& Lat() { return lat; }
   inline const T& Lon() const { return lon; }
   inline       T& Lon() { return lon; }
   friend std::ostream& operator << (std::ostream& o, Point a) {
      //o << "(" << a.lon << ", " << a.lat<< ")"; 
      o.setf(std::ios::fixed);
      o << std::left << std::setprecision(4) << a.lon << " " << a.lat;
      return o;
   }
};
#endif

template <class T>
class Path {
public:
   Path(const T lo1in=NaN, const T la1in=NaN, 
        const T lo2in=NaN, const T la2in=NaN,
		  const T dis=NaN, const T a1=NaN, const T a2=NaN)
      : long1(lo1in), lati1(la1in), long2(lo2in), lati2(la2in)
      , dist(dis), distf(dis), alpha1(a1), alpha2(a2)
		, pi(3.1415926535898), pio180(0.0174532925199433)
		, Ra(6378.137), Rb(6356.752314245), Raaobbmo(0.0067394967423330954)
		, f(0.00335281066474748), Rboa(0.99664718933525252), toler(1.0e-8) {
		if(long1 < 0.) long1 += 360.;
		if(long2 < 0.) long2 += 360.;
	}

   Path(const Point<T>& p1, const Point<T>& p2)
		: Path( p1.Lon(), p1.Lat(), p2.Lon(), p2.Lat() ) {}

   Path(const Point<T>& p1, const float dis, const float alp1)
		: Path( p1.Lon(), p1.Lat(), NaN, NaN, dis, alp1 ) {}

   ~Path(){}

	void reset(const T lo1in=NaN, const T la1in=NaN,
				  const T lo2in=NaN, const T la2in=NaN,
				  const T dis=NaN, const T a1=NaN, const T a2=NaN) {
		long1 = lo1in; lati1 = la1in; long2 = lo2in; lati2 = la2in;
      dist = distf = dis; alpha1 = a1; alpha2 = a2;
		if(long1 < 0.) long1 += 360.;
		if(long2 < 0.) long2 += 360.;
	}

	void reset(const Point<T>& p1, const Point<T>& p2) {
		reset( p1.Lon(), p1.Lat(), p2.Lon(), p2.Lat() );
	}

	void reset(const Point<T>& p1, const float dis, const float alp1) {
		reset( p1.Lon(), p1.Lat(), NaN, NaN, dis, alp1 );
	}

	const Point<T> P1() const { return Point<T>(long1, lati1); }
	const Point<T> P2() const { 
		if( long2==NaN || lati2==NaN ) calc_P2();
		return Point<T>(long2, lati2); 
	}
	const Point<T> PC() const {
		if( long2==NaN || lati2==NaN ) calc_P2();
		return Point<T>(0.5*(long1+long2), 0.5*(lati1+lati2)); 
	}

	friend std::ostream& operator << ( std::ostream& o, class Path p ) {
		o<<p.long1<<" "<<p.lati1<<"  "<<p.long2<<" "<<p.lati2<<"  "<<p.dist<<" "<<p.alpha1<<" "<<p.alpha2;
		return o;
	}

	inline const T& DistF() const {
      if( distf == NaN ) calc_dist_fast();
		return distf;
	}

   inline const T& Dist() const {
      if( dist == NaN ) calc_dist();
      return dist;
   }

   inline const T& Azi1() const {
      if( alpha1 == NaN ) calc_dist();
      return alpha1;
   }

   inline const T& Azi2() const {
      if( alpha2 == NaN )
			if( dist!=NaN && alpha1!=NaN )
				calc_P2();
			else
				calc_dist();
      return alpha2;
   }


protected:
/*
	static constexpr float NaN = -12345.;
	static constexpr double pi = 3.1415926535898;
	static constexpr double pio180 = 0.0174532925199433;
	static constexpr double Ra = 6378.137;
	static constexpr double Rb = 6356.752314245;
	static constexpr double Raaobbmo = 0.0067394967566;
	static constexpr double f = 0.00335281066474748;
	static constexpr double Rboa = 0.99664718933525252;
*/

private:
   T lati1, long1;
	mutable T lati2, long2;
   mutable T distf, dist, alpha1, alpha2;
	static const int NaN = -12345;
	const double pi;
	const double pio180;
	const double Ra;
	const double Rb;
	const double Raaobbmo;
	const double f;
	const double Rboa;
	const double toler;
	// toler ~ dis / R: 1.6e-10 toler ~ 1 mm
	mutable int iter = NaN;	// for debug

	void calc_dist_fast() const {
		if( lati1==NaN || long1==NaN || lati2==NaN || long2==NaN )
			throw std::runtime_error("Error(calc_dist_fast): empty location(s)!");

		const T dlong = fabs(long1 - long2);
		if( lati1==lati2 && (dlong==0. || dlong==360.) ) {
			distf = 0.; return;
		}

		const T ff = (lati1 + lati2) * pio180 * 0.5;
		const T g = (lati1 - lati2) * pio180 * 0.5;
		const T lambda = (long1 - long2) * pio180 * 0.5;
		const T sin_g = sin( g), cos_g = cos( g);
		const T sin_f = sin(ff), cos_f = cos(ff);
		const T sin_lambda = sin( lambda), cos_lambda = cos( lambda);
		const T sin_g2 = sin_g * sin_g, cos_g2 = cos_g * cos_g;
		const T sin_f2 = sin_f * sin_f, cos_f2 = cos_f * cos_f;
		const T sin_lambda2 = sin_lambda * sin_lambda;
		const T cos_lambda2 = cos_lambda * cos_lambda;
		const T s = sin_g2 * cos_lambda2 + cos_f2 * sin_lambda2;
		const T c = cos_g2 * cos_lambda2 + sin_f2 * sin_lambda2;
		const T omega = atan( sqrt( s / c));
		const T r = sqrt( s * c) / omega;
		const T d = 2. * omega;
		const T h1 = (3. * r - 1.) / (2. * c);
		const T h2 = (3. * r + 1.) / (2. * s);
		distf = d * (1. + f * (h1*sin_f2*cos_g2 - h2*cos_f2*sin_g2) ) * Ra;
	}

	inline static const T cvt_lat( const T b, const T lat) {
		return ( atan2( b * sin(lat), cos(lat)) );
	}

	inline static const T big_a_poly( const T usquared) {
		return( 1. + usquared *
				(4096 + usquared * (-768 + usquared * (320 - 175 * usquared))) / 16384.);
	}

	inline static const T big_b_poly( const T usquared) {
		return( usquared *
				(256. + usquared * (-128. + usquared * (74 - 47 * usquared))) / 1024.);
	}

	void calc_dist() const {
		if( lati1==NaN || long1==NaN || lati2==NaN || long2==NaN )
			throw std::runtime_error("Error(calc_dist): empty location(s)!");

		T dlt_long = fabs(long1-long2);
		if( dlt_long > 180. ) dlt_long = 360. - dlt_long;
		dlt_long *= pio180;

		// same point
		if( lati1==lati2 && dlt_long==0. ) {
			dist = distf = 0.;
			alpha1 = alpha2 = 999.;
			return ;
		}

		T U1 = cvt_lat( Rboa, lati1*pio180 );
		T U2 = cvt_lat( Rboa, lati2*pio180 );
		T sinU1 = sin(U1), cosU1 = cos(U1);
		T sinU2 = sin(U2), cosU2 = cos(U2);

		T numbda = dlt_long, numbda_old;
		int itermax = 20; //bool flag = false;
		iter = 0;
		T sigma, sin_sigma, cos_sigma;
		T cos_s_alp, cos_2sigmam;
		do {
			iter++;
			numbda_old = numbda;
			T sinnumbda = sin(numbda), cosnumbda = cos(numbda);
			T ftmp1 = cosU2*sinnumbda, ftmp2 = cosU1*sinU2-sinU1*cosU2*cosnumbda;
			T ftmp3 = sinU1*sinU2, ftmp4 = cosU1*cosU2;
			sin_sigma = sqrt( ftmp1*ftmp1 + ftmp2*ftmp2 );
			cos_sigma = ftmp3 + ftmp4*cosnumbda;
			sigma = atan2(sin_sigma,cos_sigma);
			if(sigma==0) sigma = 1.e-10;//0.00000000001;
			T sin_alp = ftmp4*sinnumbda/sin_sigma;
			cos_s_alp = 1 - sin_alp*sin_alp;
			if(cos_s_alp==0) cos_s_alp = 1.e-10;//0.00000000001;
			cos_2sigmam = cos_sigma - 2*ftmp3/cos_s_alp;
			T C = f*0.0625*cos_s_alp*(4 + f*(4 - 3*cos_s_alp));
			numbda = dlt_long + (1-C)*f*sin_alp*(sigma + C*sin_sigma*(cos_2sigmam + C*cos_sigma*(-1 +2*cos_2sigmam*cos_2sigmam)));
		} while( fabs(numbda - numbda_old)>toler && iter<itermax);

		if( iter > itermax ) {	// anti podal
			calc_dist_fast();
			dist = distf;
			alpha1 = alpha2 = 999.;
		} else {
			T mius = cos_s_alp*Raaobbmo;
			T A = big_a_poly(mius);
			T B = big_b_poly(mius);
			T d_sigma = B*sin_sigma*(cos_2sigmam +0.25*B*(cos_sigma*(-1 + 2*cos_2sigmam*cos_2sigmam)
							-B/6*cos_2sigmam*(-3+4*sin_sigma*sin_sigma)*(-3+4*cos_2sigmam*cos_2sigmam) ));
			distf = dist = Rb * A *(sigma - d_sigma);
		}

		alpha1=atan2(cosU2*sin(numbda), cosU1*sinU2 - sinU1*cosU2*cos(numbda))/pio180;
		alpha2=atan2(cosU1*sin(numbda), -sinU1*cosU2 + cosU1*sinU2*cos(numbda))/pio180;
		if( fabs(long2-long1)>180 ) { alpha1 = 360 - alpha1; alpha2 = 360 - alpha2; }
		if( long2 < long1 ) { alpha1 = 360 - alpha1; alpha2 = 360 - alpha2; }
	}

	/*
	void calc_azimuth() const {
		if( lati1==NaN || long1==NaN || lati2==NaN || long2==NaN )
			throw std::runtime_error("Error(calc_azimuth): empty location(s)!");


		T dlt_long = fabs(long1-long2);

		if(lati1==-lati2 && dlt_long==180){
			alpha1 = alpha2 = 999.;
			return;
		}

		if( dlt_long > 180. ) dlt_long = 360. - dlt_long;

		T U1 = atan(Rboa*tan(lati1*pio180));
		T U2 = atan(Rboa*tan(lati2*pio180));
		T sinU1 = sin(U1), sinU2 = sin(U2);
		T cosU1 = cos(U1), cosU2 = cos(U2);
		dlt_long = dlt_long*pio180;

		alpha1=atan2(cosU2*sin(dlt_long), cosU1*sinU2 - sinU1*cosU2*cos(dlt_long))/pio180;
		alpha2=atan2(cosU1*sin(dlt_long), -sinU1*cosU2 + cosU1*sinU2*cos(dlt_long))/pio180;
		if( fabs(long2-long1)>180 ) { alpha1 = 360 - alpha1; alpha2 = 360 - alpha2; }
		if( long2 < long1 ) { alpha1 = 360 - alpha1; alpha2 = 360 - alpha2; }

	} 
	*/

	void calc_P2() const {
		if( lati1==NaN || long1==NaN || dist==NaN || alpha1==NaN )
			throw std::runtime_error("Error(calc_P2): record incomplete");
		const T Rboa = 1. - f;       /* if a = 1 */
		/* u1, u2 are "reduced latitudes" */
		const T latio1=lati1*pio180, alphao1=alpha1*pio180;
		const T u1 = cvt_lat( Rboa, latio1);
		const T cos_u1 = cos(u1), sin_u1 = sin(u1);
		const T sin_a1 = sin(alphao1), cos_a1 = cos(alphao1);
		const T sigma1 = atan2( sin_u1, cos_a1 * cos_u1 );
		const T sin_alpha = cos_u1 * sin_a1;
		const T cos2_alpha = 1. - sin_alpha * sin_alpha;
		const T u_squared = cos2_alpha * ( 1./(Rboa*Rboa) - 1. );
		const T big_a = big_a_poly(u_squared);
		const T big_b = big_b_poly( u_squared);

		T sigma0 = dist / (Rb * big_a);
		T sigma = sigma0, cos_2sigma_m, sin_sigma, cos_sigma;
		T sigma_old = sigma + toler*2;
		do {
			cos_2sigma_m = cos(2. * sigma1 + sigma);
			T cos_2sigma_m2 = cos_2sigma_m * cos_2sigma_m;
			sin_sigma = sin(sigma);
			T delta_sigma = big_b * sin_sigma * ( cos_2sigma_m + (0.25*big_b) * ( cos(sigma)*(2.*cos_2sigma_m2 - 1.) - 
								 (big_b/6.)*cos_2sigma_m*(4.*sin_sigma*sin_sigma-3.)*(4.*cos_2sigma_m2-3.) ) );
			sigma_old = sigma;
			sigma = sigma0 + delta_sigma;
		} while ( (sigma - sigma_old) > toler );

		sin_sigma = sin(sigma);
		cos_sigma = cos(sigma);
		cos_2sigma_m = cos(2. * sigma1 + sigma);
		const T tval = sin_u1 * sin_sigma - cos_u1 * cos_sigma * cos_a1;
		lati2 = atan2(sin_u1 * cos_sigma + cos_u1 * sin_sigma * cos_a1,
				  Rboa * sqrt(sin_alpha * sin_alpha + tval * tval)) / pio180;
		T lambda = atan2(sin_sigma * sin_a1, cos_u1*cos_sigma - sin_u1*sin_sigma*cos_a1);
		// T big_c = f / (16. * cos2_alpha * (4. + f * (4. - 3 * cos2_alpha)));
		T big_c = (f / 16) * cos2_alpha * (4. + f * (4. - 3 * cos2_alpha));
		T big_l = lambda - (1. - big_c) * f * sin_alpha *
					 (sigma + big_c * sin_sigma * (cos_2sigma_m + big_c * cos_sigma *
					 (2. * cos_2sigma_m * cos_2sigma_m - 1.)));
		long2 = long1 + big_l/pio180;
		alpha2 = atan2(sin_alpha, cos_u1*cos_sigma*cos_a1 - sin_u1*sin_sigma)/pio180;
		if( alpha2 < 0. ) alpha2 += 360.;
	}

};

/*
	template <class T>
	Point<T> PgeoToPxy( const Point<T> Pgeo, const Point<T> Pcenter ) {
	auto path = Path<T>(Pcenter, Pgeo);
	T dis = path.Dist(), azi = path.Azi1()*M_PI/180.;
	return Point<T>( dis*sin(azi), dis*cos(azi) );
	}


	template <class T>
	bool Check_azi_cov(T *lon, T *lat, int nsta, T mdis, T *flag) {
	int i, j, jj, k, ndata;
	T azi[300], azimin, azitmp, dist;
	T mazi=110.;

	for(i=0;i<nsta;i++){
	if(flag[i]==-1) continue;
	ndata=0;
	for(j=0;j<nsta;j++){
	if(j==i || flag[j]==-1) continue;
	calc_dist(lat[i], lon[i], lat[j], lon[j], &dist);
	if(dist>mdis) continue;
	calc_azimuth(lat[i], lon[i], lat[j], lon[j], &azi[ndata]);
	ndata++;
	}
//cout<<"Check azi cov: "<<ndata<<" nearby stations"<<endl;
for(j=0;j<ndata;j++) {
azimin=azi[j]; jj=j;
for(k=j+1;k<ndata;k++)
if(azimin>azi[k]) {
azimin=azi[k];
jj=k;
}
if(jj==j) continue;
azitmp=azi[j];
azi[j]=azi[jj];
azi[jj]=azitmp;
}
if(azi[0]+360.-azi[ndata-1]>mazi) {
flag[i]=-2;
continue;
}
			for(j=1;j<ndata;j++)
				if(azi[j]-azi[j-1]>mazi) {
					flag[i]=-2;
					break;
				}
		}
		for(i=0;i<nsta;i++) if(flag[i]==-2) flag[i]=-1;
		return 1;
	}
*/

#endif
