#include "SacRec.h"
#include <iostream>
#include <fstream>

class Parabola;
class Point{
public:
	Point(const float xin = NaN, const float yin = NaN)
		: x(xin), y(yin) {}

	friend std::ostream& operator<<(std::ostream& o, const Point& p) {
		o<<p.x<<" "<<p.y;
		return o;
	}

	friend Parabola;

protected:
	static constexpr float NaN = -12345.;

private:
	float x, y;
};

class Parabola {
public:
	Parabola( const Point& P1in, const Point& P2in, const Point& P3in )
		: P1(P1in), P2(P2in), P3(P3in) {}

	void Solve(){
		float x1 = P1.x, y1 = P1.y;
		float x2 = P2.x, y2 = P2.y;
		float x3 = P3.x, y3 = P3.y;
		float xs1 = x1*x1, xs2 = x2*x2, xs3 = x3*x3;
		float denom = (x1-x2) * (x1-x3) * (x2-x3);
		a = (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)) / denom;
		b = (xs3*(y1-y2) + xs2*(y3-y1) + xs1*(y2-y3)) / denom;
		c = (x2*x3*(x2-x3)*y1 + x3*x1*(x3-x1)*y2 + x1*x2*(x1-x2)*y3) / denom;
	}

	Point Vertex() {
		if(a==NaN || b==NaN || c==NaN)
			Solve();
		if( PV.x==NaN || PV.y==NaN ) {
			PV.x = - b / (2.*a);
			PV.y = c - a*PV.x*PV.x;
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


int main( int argc, char* argv[] ) {
	if( argc != 4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sac1] [sac2] [fout]"<<std::endl;
		return -1;
	}

	// correlate
	SacRec sac1, sac2, saccor;
	sac1.Load(argv[1]);
	sac2.Load(argv[2]);
	sac1.Correlate(sac2, saccor);

	// file for output 
	std::ofstream fout(argv[3]);
	if( ! fout ) {
		std::cerr<<"Error(main): cannot write to file "<<argv[3]<<std::endl;
		return -2;
	}

	// find T of maximum correlation at each period
	float perfactor = 1.05;
	float shdb = saccor.shd.b, delta = saccor.shd.delta;
	for( float per=2.; per<50.; per*=perfactor ) {
		// gauss filt
		SacRec sacflt;
		saccor.GaussianFilt( 1./per, 0.02, sacflt );
		// find maximum
		int imin, imax;
		sacflt.MinMax(imin, imax, -per*0.51, per*0.51);
		if( imax<=0 || imax> sacflt.shd.npts-2 ) continue;
		// the max, left, and right points
		float *sigsac = sacflt.sig.get();
		Point p1( shdb + (imax-1)*delta, sigsac[imax-1] );
		Point p2( shdb + imax*delta, sigsac[imax] );
		Point p3( shdb + (imax+1)*delta, sigsac[imax+1] );
		Parabola pr(p1, p2, p3);
		fout<<per<<" "<<pr.Vertex()<<std::endl;
	}
	
}
