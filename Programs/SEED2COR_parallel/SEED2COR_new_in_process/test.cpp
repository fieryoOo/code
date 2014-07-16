#include <iostream>

struct X {
	float datax;
};

struct Y {
	float datay;
};
struct XYZ : public X, public Y {
	float dataz;
	XYZ() {}
	XYZ( const X& x ) : X(x) {}
	XYZ( const Y& y ) : Y(y) {}
};

int main() {
	X x;
	x.datax = 123.;
	Y y;
	y.datay = 456.;
	XYZ xyz;
	xyz = x;
	std::cerr<<xyz.datax<<std::endl;
	xyz = y;
	std::cerr<<xyz.datax<<std::endl;
}
