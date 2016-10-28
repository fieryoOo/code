#ifndef PARABOLA_H
#define PARABOLA_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdexcept>

#ifndef FuncName
#define FuncName __FUNCTION__
#endif

#ifndef POINTC
#define POINTC
/* ----- single point data structures ----- */
//template <class T> class Curve;
//class KDeriv;
//class Parabola;
class PointC {
public:
	PointC( float xin=NaN, float yin=NaN, float zin=1. ) 
		: x(xin), y(yin), z(zin) {
//std::cerr<<" PointC::PointC: "<<zin<<std::endl;
	}

	PointC( const std::string& input ) {
		int nrd = sscanf(input.c_str(), "%f %f %f", &x, &y, &z);
		if( nrd < 2 )
			throw std::runtime_error( std::string(FuncName) + ": format error in string "+input );
	}

	bool isValid() const { return (x!=NaN && y!=NaN); }

	// Norm
	float abs() const { return sqrt(x*x + y*y); }

	// unary negation 
	PointC operator-() const { return PointC( -x, -y, -z ); }

	// addition 
	PointC& operator+=( const PointC& p2 ) {
		x += p2.x; y += p2.y; z += p2.z;
		return *this; 
	}
	friend PointC operator+( const PointC& p1, const PointC& p2 ) {
		PointC pres = p1;
		pres += p2;
		return pres;
	}

	// subtraction
   PointC& operator-=( const PointC& p2 ) {
		(*this) += -p2;
      return *this;
   }
   friend PointC operator-( const PointC& p1, const PointC& p2 ) {
		PointC pres = p1;
      pres -= p2;
      return pres;
   }

	// multiplication (*float)
	PointC& operator*=( const float mul ) { 
		x *= mul; y *= mul; z *= mul;
		return *this; 
	}
	friend PointC operator*( const PointC& p1, float mul ) {
		PointC pres = p1;
		pres *= mul;
		return pres;
	}
	friend PointC operator*( float mul, const PointC& p1 ) {
		PointC pres = p1;
		pres *= mul;
		return pres;
	}

	friend bool operator< ( const PointC& p1, const PointC& p2 ) {
		return (p1.x < p2.x);
	}

	// point - line distance
	friend float PLDist( const PointC& p, const PointC& pl1, const PointC& pl2 ) {
		const PointC pt1 = p - pl1, pt2 = pl2 - pl1;
		return fabs( pt1.x*pt2.y - pt2.x*pt1.y ) / pt2.abs();
	}

	friend std::ostream& operator<< ( std::ostream& o, const PointC& pt ) {
		o << pt.x << " " << pt.y << " " << pt.z;
		return o;
	}

	//friend Curve<PointC>;
	//friend KDeriv;
	//friend Parabola;
	//friend Curve<PointC> operator-(const Curve<PointC>& c1, const Curve<PointC>& c2);

	static constexpr float NaN = -123456.;
	static constexpr float twopi = M_PI * 2.;

//protected:
	float x = NaN, y = NaN, z = 1.;
};
#endif	// POINTC


// a parabola of the standard form y = ax**2 + bx + c
class Parabola {
public:
	Parabola() {}
	Parabola( const PointC& P1in, const PointC& P2in, const PointC& P3in )
		: P1(P1in), P2(P2in), P3(P3in) {}

	void Solve() const {
		if( ! (P1.isValid() && P2.isValid() && P3.isValid()) ) {
			std::stringstream ss; ss<<FuncName<<": unfilled point(s). "<<P1<<"   "<<P2<<"   "<<P3;
			throw std::runtime_error( ss.str() );
		}
		long double x1 = P1.x, y1 = P1.y;
		long double x2 = P2.x, y2 = P2.y;
		long double x3 = P3.x, y3 = P3.y;
		long double xs1 = x1*x1, xs2 = x2*x2, xs3 = x3*x3;
		long double denom = (x1-x2) * (x1-x3) * (x2-x3);
		a = (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)) / denom;
		b = (xs3*(y1-y2) + xs2*(y3-y1) + xs1*(y2-y3)) / denom;
		c = (x2*x3*(x2-x3)*y1 + x3*x1*(x3-x1)*y2 + x1*x2*(x1-x2)*y3) / denom;
		//std::cerr<<"results: "<<a<<" "<<b<<" "<<c<<std::endl;
	}

	PointC Vertex() const {
		if(a==NaN || b==NaN || c==NaN) Solve();
		if( PV.x==NaN || PV.y==NaN ) {
			PV.x = - b / (2.*a);
			PV.y = c - a*PV.x*PV.x;
		}
		return PV;
	}

	long double A() const { 
		if(a==NaN) Solve(); 
		return a;
	}
	long double B() const { 
		if(b==NaN) Solve(); 
		return b;
	}
	long double C() const { 
		if(c==NaN) Solve(); 
		return c;
	}

	long double Slope( const long double x ) const {
		if(a==NaN || b==NaN || c==NaN) Solve();
		return ( 2*a*x + b );
	}

	long double operator[]( const long double x ) const {
		return ( a*x*x + b*x + c );
	}
	long double operator()( const long double x ) const {
		if(a==NaN || b==NaN || c==NaN) Solve();
		return ( a*x*x + b*x + c );
	}

	// compute least-fit parabola to dataV, a (weighted) rms is also computed
	void Fit(typename std::vector<PointC>::const_iterator iterl, 
				typename std::vector<PointC>::const_iterator iteru, float& rmsw) {
		// allocate working space
		int dim = 3, ndat = iteru - iterl;	//dataV.size();
		if(ndat < 3) {
			printf("no enough data points for fitting");
			return ;
		} else if(ndat == 3) {
			P1 = *iterl; P2 = *(iterl+1); P3 = *(iterl+2);
			Solve(); return;
		}
		double A[dim][ndat], AC[dim];
		double **AA = (double **) malloc ( dim * sizeof(double *) );
		double **AAI = (double **) malloc ( dim * sizeof(double *) );
		AA[0] = (double*) malloc ( dim * dim * sizeof(double) );
		AAI[0] = (double*) malloc ( dim * dim * sizeof(double) );
		for(int i=1;i<dim;i++) {
			AA[i] = AA[i-1] + dim; 
			AAI[i] = AAI[i-1] + dim;
		}

		// fill matrix A
	   for(int i=0;i<ndat;i++) {
			const auto& p = *(iterl+i);	//dataV[i];
			float x = p.x, w = p.z;
	      A[0][i] = x*x*w;
	      A[1][i] = x*w;
	      A[2][i] = w;
	   }

		// compute ATA
		for(int i=0;i<dim;i++)
			for(int j=i;j<dim;j++) {
				AA[i][j] = 0.;
				for(int ii=0;ii<ndat;ii++) AA[i][j] += A[i][ii]*A[j][ii];
			}
		for(int i=1;i<dim;i++) for(int j=0;j<i;j++) AA[i][j] = AA[j][i];
		// and its inverse
		Inverse( AA, dim, AAI );

		// compute AC
		for(int i=0;i<dim;i++) {
			AC[i]=0;
			for(int ii=0;ii<ndat;ii++) {
				const auto& p = *(iterl+ii);	//dataV[ii];
				float y = p.y, w = p.z;
				AC[i] += A[i][ii] * y * w;
			}
		}

		// fit results
		double coef[dim];
	   for(int i=0;i<dim;i++) {
		   coef[i] = 0;
			for(int j=0;j<dim;j++) {
				coef[i] += AAI[i][j]*AC[j];
			}
	   }
		a = coef[0];
		b = coef[1];
		c = coef[2];

		// clean up
		free(AA[0]); free(AAI[0]);
		free(AA); free(AAI);

		// errors
		rmsw = 0; //std = 0;
		double v1=0, v2=0;
		for(int i=0;i<ndat;i++) {
			const auto& p = *(iterl+i);	//dataV[i];
			float x = p.x, y = p.y, w = p.z;
			float diff = y - (a * x*x + b * x + c);
			//std = std + diff*diff;
			rmsw += diff*diff * w;
			v1 += w;
			v2 += w * w;
		}
		//std = sqrt(std/(ndat-3.));
		float denom = v1*v1-3.*v2;
		if( denom <= 0. ) {
			rmsw = NaN;
		} else {
			rmsw = sqrt(rmsw*v1/denom);
		}
		//cout<<*A0<<" + "<<*A1<<" * sin( theta + "<<*phi<<" )"<<endl;
	}
	void Fit( const std::vector<PointC>& dataV, float& rmsw ) {
		Fit( dataV.begin(), dataV.end(), rmsw );
	}

	friend std::ostream& operator<<( std::ostream& o, const Parabola& parab) {
		o<<parab.a<<" "<<parab.b<<" "<<parab.c;
		return o;
	}

protected:
	static const int NMAX = 10;
	static constexpr long double NaN = PointC::NaN;

private:
	PointC P1, P2, P3;
	mutable PointC PV;
	mutable long double a=NaN, b=NaN, c=NaN;

	void arg(double *a, double *b, int *n,int x,int y) {
		int k,l,i,j;
		for(i=0,k=0;i<*n;i++,k++) {
			for(j=0,l=0;j<*n;j++,l++) {
				if(i==x)
					i++;
				if(j==y)
					j++;
				*(b+NMAX*k+l)=*(a+NMAX*i+j);
			}
		}
		*n=*n-1;
	}

	double det(double *p,int *n) {
		int i,j,m;
		double d[NMAX][NMAX], sum=0;
		m=*n;
		if(*n==2)
			return(*p**(p+NMAX+1)-*(p+1)**(p+NMAX));
		for(i=0,j=0;j<m;j++) {
			*n=m;
			arg(p,&d[0][0],n,i,j);
			sum=sum+*(p+NMAX*i+j)*pow(-1,(i+j))*det(&d[0][0],n);
		}

		return(sum);
	}

	int Inverse( double **datin, int n, double **datout )	{
		//void arg(int *,int *, int *,int ,int );
		//int det(int *,int *);
		int i,j,m;
		double a[NMAX][NMAX],b[NMAX][NMAX],c[NMAX][NMAX],d;
		//clrscr();
		for(i=0;i<n;i++) for(j=0;j<n;j++) a[i][j] = datin[i][j];
		if(n==2) {
			c[0][0]=a[1][1];
			c[1][1]=a[0][0];
			c[0][1]=-a[0][1];
			c[1][0]=-a[1][0];
			d=a[0][0]*a[1][1]-a[0][1]*a[1][0];
			//printf("Determinant: %lf\n",d);
			if(d==0) {
				//getch();
				//std::cin.get();
				return 0; //exit((int)d-'0');
			}

			for(i=0;i<n;i++) {
				printf("\n");
				for(j=0;j<n;j++)
					printf(" %f",c[i][j]/(float)d);
			}
		}
		else {
			m=n;
			for(i=0;i<m;i++) {
				for(j=0;j<m;j++) {
					n=m;
					arg(&a[0][0],&b[0][0],&n,i,j);
					c[j][i]=pow(-1,(i+j))*det(&b[0][0],&n);
				}
			}
			n=m;
			d=det(&a[0][0],&n);
			//printf("Determinant is :%d\n",d);
			if(d==0) {
				printf("INVERSE DOES NOT EXIST");
				//getch();
				//std::cin.get();
				return 0; //exit((int)d-'0');
			}
			for(i=0;i<m;i++) {
				for(j=0;j<m;j++)
					datout[i][j] = c[i][j]/d; //printf(" %f",c[i][j]/d);
			}
		} //std::cin.get();
		return 1;
	}

};

#endif
